#!/bin/bash
#SBATCH --job-name=betageo-gen
#SBATCH --time=04:00:00
#SBATCH --mem=4G
#SBATCH --partition=batch
#SBATCH --output=logs/gen_%A_%a.out
#SBATCH --error=logs/gen_%A_%a.err
#SBATCH --array=0-127

# Split the seeds first, then set --array above to match what split_seeds.py prints:
#   python split_seeds.py --range 0-9999 --tasks 128

module load python/3.11.0s-ixrhc3q
source ~/geo_env/bin/activate
cd ~/betageometry

mkdir -p logs data

BATCH_FILE="seed_batches/batch_${SLURM_ARRAY_TASK_ID}.txt"
if [ ! -f "$BATCH_FILE" ]; then
    echo "Error: $BATCH_FILE does not exist. Did you run split_seeds.py?"
    exit 1
fi

echo "Task $SLURM_ARRAY_TASK_ID processing $BATCH_FILE"
python generate_batch.py "$BATCH_FILE" --quiet
