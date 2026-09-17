import json
import random
import torch
from datasets import Dataset
from transformers import (
    AutoTokenizer, 
    AutoModelForCausalLM, 
    Trainer, 
    TrainingArguments, 
    DataCollatorForLanguageModeling
)
from peft import LoraConfig, get_peft_model

# --- GLOBALS ---
# File paths
JSON_FILE = "data/training_data_clean.json"
OUTPUT_DIR = "./gemma-finetuned-geometry"

# Model ID
MODEL_ID = "google/gemma-3-1b-pt"

# Data hyperparameters
SEED = 42
SUBSET_SIZE = 100000
MAX_SEQ_LENGTH = 512
SEPARATOR = "\n### Solution:\n"
INSTRUCTION = "Use a new point name not listed below as the first point in the construction. Only output one construction."

# Training hyperparameters
NUM_EPOCHS = 5           
LEARNING_RATE = 2e-4
BATCH_SIZE = 2
GRAD_ACCUMULATION = 2
LORA_R = 8
LORA_ALPHA = 32
LORA_DROPOUT = 0.1

# ---LOAD DATA ---
random.seed(SEED)

with open(JSON_FILE, "r") as f:
    data = json.load(f)

# Take a random subset
subset_data = random.sample(data, min(SUBSET_SIZE, len(data)))
records = [{"input_text": d["input_text"], "output_text": d["output_text"]} for d in subset_data]
dataset = Dataset.from_list(records)

print(f"Loaded {len(dataset)} examples.")

# --- LOAD TOKENIZER ---
tokenizer = AutoTokenizer.from_pretrained(MODEL_ID)
tokenizer.pad_token = tokenizer.eos_token

# --- TOKENIZATION FUNCTION (with Instruction & Masking) ---
def tokenize_fn(examples):
    model_inputs = {"input_ids": [], "attention_mask": [], "labels": []}
    
    for inp, out in zip(examples["input_text"], examples["output_text"]):
        # Structure: <s> Input + Instruction + Separator + Output + <eos>
        full_text = f"<s>{inp} {INSTRUCTION}{SEPARATOR}{out}{tokenizer.eos_token}"
        
        # Tokenize the full sequence
        tokenized = tokenizer(
            full_text,
            truncation=True,
            max_length=MAX_SEQ_LENGTH,
            padding="max_length",
            return_tensors=None
        )
        
        input_ids = tokenized["input_ids"]
        labels = input_ids.copy()

        # We re-tokenize the PROMPT part (Input + Instruction + Separator)
        prompt_part = f"<s>{inp} {INSTRUCTION}{SEPARATOR}"
        prompt_len = len(tokenizer(prompt_part, add_special_tokens=False)["input_ids"])

        for i in range(len(labels)):
            # Mask the input prompt (which now includes the instruction) AND padding
            if i < prompt_len or labels[i] == tokenizer.pad_token_id:
                labels[i] = -100
            # Otherwise leave the label as the token ID (for the answer part)

        model_inputs["input_ids"].append(input_ids)
        model_inputs["attention_mask"].append(tokenized["attention_mask"])
        model_inputs["labels"].append(labels)
        
    return model_inputs

# Apply tokenization
tokenized_dataset = dataset.map(tokenize_fn, batched=True, remove_columns=dataset.column_names)

# --- LOAD MODEL & APPLY LoRA ---
model = AutoModelForCausalLM.from_pretrained(
    MODEL_ID,
    dtype=torch.bfloat16,
    device_map="auto"
)

lora_config = LoraConfig(
    r=LORA_R,
    lora_alpha=LORA_ALPHA,
    target_modules=["q_proj", "v_proj"],
    lora_dropout=LORA_DROPOUT,
    bias="none",
    task_type="CAUSAL_LM"
)

model = get_peft_model(model, lora_config)
model.print_trainable_parameters()

# --- SETUP TRAINER ---
data_collator = DataCollatorForLanguageModeling(tokenizer=tokenizer, mlm=False)

training_args = TrainingArguments(
    output_dir=OUTPUT_DIR,
    per_device_train_batch_size=BATCH_SIZE,
    gradient_accumulation_steps=GRAD_ACCUMULATION,
    learning_rate=LEARNING_RATE,
    num_train_epochs=NUM_EPOCHS,
    logging_steps=5,
    save_strategy="epoch",
    save_total_limit=1,
    bf16=True,
    optim="adamw_torch",
    report_to="none"
)

trainer = Trainer(
    model=model,
    train_dataset=tokenized_dataset,
    args=training_args,
    data_collator=data_collator
)

# --- TRAIN AND SAVE MODEL ---
print("Starting fine-tuning...")
trainer.train()
trainer.save_model(OUTPUT_DIR)
print(f"Fine-tuning complete. Model saved to {OUTPUT_DIR}")

# --- INFERENCE TEST ---
print("\n--- Running Inference Test ---")
model.eval()

# Example geometry problem
raw_prompt = "midp X9 X1 X4; midp X8 X1 X2; col X5 X1 X2; perp X4 X5 X1 X2; midp X7 X4 X2; ? cong X9 X1 X8 X7"

# Construct the prompt exactly how it was trained: Input + Instruction + Separator
formatted_prompt = f"{raw_prompt} {INSTRUCTION}{SEPARATOR}"

inputs = tokenizer(formatted_prompt, return_tensors="pt").to(model.device)

with torch.inference_mode():
    outputs = model.generate(
        **inputs,
        max_new_tokens=50,
        do_sample=False,
        eos_token_id=tokenizer.eos_token_id,
        pad_token_id=tokenizer.eos_token_id
    )

generated_text = tokenizer.decode(outputs[0][inputs["input_ids"].shape[-1]:], skip_special_tokens=True)

print(f"Input Prompt: {raw_prompt}")
print(f"Generated Answer: {generated_text.strip()}")