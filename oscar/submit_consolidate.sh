#!/bin/bash
#SBATCH --job-name=geo_merge
#SBATCH --time=01:00:00
#SBATCH --mem=16G                 # More memory to handle large JSON merge
#SBATCH --partition=batch
#SBATCH --output=logs/merge_%j.out
#SBATCH --error=logs/merge_%j.err

# Load Environment
module load python/3.9.16s-x3wdtvt
source ~/geo_env/bin/activate
cd ~/betageometry

echo "Starting data consolidation..."
python training_file_generation.py
echo "Consolidation complete."