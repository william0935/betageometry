"""LoRA fine-tuning of Gemma on generated construction data.

Ported from BetaGeometry's `finetuning_model.py`, with three corrections. Each of them is
silent -- the script trains and reports a falling loss either way -- so they are called
out here rather than just fixed:

1. The prompt was hand-built as `f"<s>{inp} ..."`. `<s>` is Llama/Mistral's BOS marker,
   not Gemma's, so it tokenised as literal text; and because `tokenizer(...)` also adds
   its own BOS by default, every example began with a stray token pair. The prompt is now
   built by `gemma.build_prompt`, which inference uses too, so the two cannot drift.

2. Label masking used `labels[i] == tokenizer.pad_token_id` to blank the padding. With
   `pad_token = eos_token` -- set two lines earlier -- that also blanks the one real EOS
   at the end of the target, so the model is never trained to stop and generation runs to
   the token limit. Padding is now identified by the attention mask instead.

3. The prompt length used to mask the question was measured with
   `add_special_tokens=False` while the full sequence was tokenised with the default
   `True`, so the boundary was off by the BOS token and one target token was masked (or
   one prompt token supervised). Both are now tokenised the same way.
"""

import argparse
import json
import random

DEFAULT_MODEL_ID = "google/gemma-3-1b-pt"


def build_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--data", default="data/training_data.json")
    p.add_argument("--output-dir", default="gemma-finetuned-geometry")
    p.add_argument("--model-id", default=DEFAULT_MODEL_ID)
    p.add_argument("--subset-size", type=int, default=0,
                   help="train on a random subset of this many examples (0 = all)")
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--max-seq-length", type=int, default=512)
    p.add_argument("--epochs", type=float, default=5)
    p.add_argument("--learning-rate", type=float, default=2e-4)
    p.add_argument("--batch-size", type=int, default=2)
    p.add_argument("--grad-accumulation", type=int, default=2)
    p.add_argument("--lora-r", type=int, default=8)
    p.add_argument("--lora-alpha", type=int, default=32)
    p.add_argument("--lora-dropout", type=float, default=0.1)
    return p.parse_args()


def load_examples(path: str, subset_size: int, seed: int):
    with open(path, "r", encoding="utf-8") as f:
        data = json.load(f)
    if not data:
        raise SystemExit(f"{path} has no examples; run data_generation.py first")
    if subset_size and subset_size < len(data):
        data = random.Random(seed).sample(data, subset_size)
    return [{"input_text": d["input_text"], "output_text": d["output_text"]} for d in data]


def main():
    args = build_args()

    import torch
    from datasets import Dataset
    from peft import LoraConfig, get_peft_model
    from transformers import (
        AutoModelForCausalLM,
        AutoTokenizer,
        DataCollatorForLanguageModeling,
        Trainer,
        TrainingArguments,
    )

    from gemma import build_example, build_prompt

    random.seed(args.seed)
    records = load_examples(args.data, args.subset_size, args.seed)
    dataset = Dataset.from_list(records)
    print(f"Loaded {len(dataset)} examples from {args.data}.")

    tokenizer = AutoTokenizer.from_pretrained(args.model_id)
    if tokenizer.pad_token is None:
        tokenizer.pad_token = tokenizer.eos_token

    def tokenize_fn(examples):
        out = {"input_ids": [], "attention_mask": [], "labels": []}
        for inp, target in zip(examples["input_text"], examples["output_text"]):
            full = build_example(inp, target, tokenizer.eos_token)
            prompt = build_prompt(inp)

            encoded = tokenizer(full, truncation=True, max_length=args.max_seq_length,
                                padding="max_length")
            # Same tokenizer settings as `full`, so the boundary lines up exactly.
            prompt_len = len(tokenizer(prompt, truncation=True,
                                       max_length=args.max_seq_length)["input_ids"])

            input_ids = encoded["input_ids"]
            attention = encoded["attention_mask"]
            labels = list(input_ids)
            for i in range(len(labels)):
                # Supervise the construction only: mask the question, and mask padding
                # via the attention mask so the target's own EOS survives and the model
                # learns where to stop.
                if i < prompt_len or attention[i] == 0:
                    labels[i] = -100

            out["input_ids"].append(input_ids)
            out["attention_mask"].append(attention)
            out["labels"].append(labels)
        return out

    tokenized = dataset.map(tokenize_fn, batched=True,
                            remove_columns=dataset.column_names)

    model = AutoModelForCausalLM.from_pretrained(
        args.model_id, dtype=torch.bfloat16, device_map="auto"
    )
    model = get_peft_model(model, LoraConfig(
        r=args.lora_r,
        lora_alpha=args.lora_alpha,
        target_modules=["q_proj", "v_proj"],
        lora_dropout=args.lora_dropout,
        bias="none",
        task_type="CAUSAL_LM",
    ))
    model.print_trainable_parameters()

    trainer = Trainer(
        model=model,
        train_dataset=tokenized,
        args=TrainingArguments(
            output_dir=args.output_dir,
            per_device_train_batch_size=args.batch_size,
            gradient_accumulation_steps=args.grad_accumulation,
            learning_rate=args.learning_rate,
            num_train_epochs=args.epochs,
            logging_steps=5,
            save_strategy="epoch",
            save_total_limit=1,
            bf16=True,
            optim="adamw_torch",
            report_to="none",
        ),
        data_collator=DataCollatorForLanguageModeling(tokenizer=tokenizer, mlm=False),
    )

    print("Starting fine-tuning...")
    trainer.train()
    trainer.save_model(args.output_dir)
    tokenizer.save_pretrained(args.output_dir)
    print(f"Fine-tuning complete. Adapter saved to {args.output_dir}")


if __name__ == "__main__":
    main()
