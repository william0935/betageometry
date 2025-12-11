#!/bin/bash
#SBATCH --job-name=betageometry
#SBATCH --time=10:00:00
#SBATCH --mem=4G
#SBATCH --partition=batch
#SBATCH --output=logs/batch_%a.out
#SBATCH --error=logs/batch_%a.err
#SBATCH --array=0-127                # 128 tasks (0 to 127)

# Load Environment
module load python/3.9.16s-x3wdtvt
source ~/geo_env/bin/activate
cd ~/betageometry

# Determine the input file for this specific task
# It looks for files named 'batch_0.txt', 'batch_1.txt', etc. in the 'seed_batches' folder
BATCH_FILE="seed_batches/batch_${SLURM_ARRAY_TASK_ID}.txt"

echo "Job ID: $SLURM_ARRAY_JOB_ID, Task ID: $SLURM_ARRAY_TASK_ID"
echo "Processing seeds from file: $BATCH_FILE"

if [ ! -f "$BATCH_FILE" ]; then
    echo "Error: Batch file $BATCH_FILE does not exist!"
    exit 1
fi

# Run the python script passing the BATCH FILE as the argument
python new_oscar_data_generation.py "$BATCH_FILE"