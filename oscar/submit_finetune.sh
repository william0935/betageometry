#!/bin/bash
#SBATCH --job-name=betageo-finetune
#SBATCH --time=08:00:00
#SBATCH --mem=32G
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --output=logs/ft_%j.out
#SBATCH --error=logs/ft_%j.err

module load python/3.11.0s-ixrhc3q
module load cuda
source ~/geo_env/bin/activate
cd ~/betageometry

# Needs the llm extra:  pip install -e '.[llm]'
# Gemma is gated on HuggingFace: accept the licence and `huggingface-cli login` first.
python finetune.py \
    --data data/training_data.json \
    --output-dir gemma-finetuned-geometry \
    --epochs 3
