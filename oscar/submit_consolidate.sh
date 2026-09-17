#!/bin/bash
#SBATCH --job-name=betageo-merge
#SBATCH --time=01:00:00
#SBATCH --mem=8G
#SBATCH --partition=batch
#SBATCH --output=logs/merge_%j.out
#SBATCH --error=logs/merge_%j.err

# Chain this after the generation array so it waits for it to finish:
#   sbatch --dependency=afterok:<ARRAY_JOB_ID> oscar/submit_consolidate.sh

module load python/3.11.0s-ixrhc3q
source ~/geo_env/bin/activate
cd ~/betageometry

python consolidate_data.py --data-dir data
