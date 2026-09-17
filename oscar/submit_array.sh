#!/bin/bash
#SBATCH --job-name=betageometry
#SBATCH --time=04:00:00           # 4 hours max per job
#SBATCH --mem=4G                  # 4GB RAM per job
#SBATCH --partition=batch
#SBATCH --output=logs/gen_%A_%a.out   # Log files
#SBATCH --error=logs/gen_%A_%a.err
#SBATCH --array=0-99              # Array range: 0 to 99 (100 total jobs)

# --- CONFIGURATION ---
BATCH_SIZE=100
START_OFFSET=90000
# ---------------------

# Calculate seed range for this specific task
START_SEED=$((START_OFFSET + SLURM_ARRAY_TASK_ID * BATCH_SIZE))
END_SEED=$((START_SEED + BATCH_SIZE))

# Load Environment
module load python/3.9.16s-x3wdtvt
source ~/geo_env/bin/activate
cd ~/betageometry

# Create directories
mkdir -p logs
mkdir -p data

echo "Task ID: $SLURM_ARRAY_TASK_ID"
echo "Processing seeds: $START_SEED to $END_SEED"

# Run the python script
python oscar_data_generation.py $START_SEED $END_SEED