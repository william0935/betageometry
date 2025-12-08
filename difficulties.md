## Rough attempt to finetune gemma

import os
import optax
import treescope
from kauldron import kd
from gemma import gm
from datasets import load_dataset

#initiallization
tokenizer = gm.text.Gemma3Tokenizer()
model = gm.nn.Gemma3_4B(tokens="batch.input",)
loss = kd.losses.SoftmaxCrossEntropyWithIntLabels(
logits="preds.logits",
labels="batch.target",
mask="batch.loss_mask",
)

#dataset -- this is the sketchiest part. Not sure how to correctly load/format our file
ds = load_dataset("json", data_files="training_data.jsonl")
def format_example(example):
if "input_text" in example and "output_text" in example:
return {
"text": f"<start_of_turn>user\n{example['input_text']}<end_of_turn>\n"
f"<start_of_turn>assistant\n{example['output_text']}<end_of_turn>"
}
return example

dataset = ds.map(format_example)

#train
trainer = kd.train.Trainer(
seed=42, # The seed of enlightenment
workdir='/tmp/ckpts', # TODO(epot): Make the workdir optional by default # Dataset
train_ds=ds, # Model
model=model,
init_transform=gm.ckpts.LoadCheckpoint( # Load the weights from the pretrained checkpoint
path=gm.ckpts.CheckpointPath.GEMMA3_4B_IT,
), # Training parameters
num_train_steps=300,
train_losses={"loss": loss},
optimizer=optax.adafactor(learning_rate=1e-3),
)

trainer.train()
