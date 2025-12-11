#!/bin/bash
#SBATCH --job-name=betageometry_cleanup
#SBATCH --time=01:00:00             # 1 hour is plenty for file scanning
#SBATCH --mem=4G                    # Low memory requirement
#SBATCH --partition=batch
#SBATCH --output=logs/clean_%j.out
#SBATCH --error=logs/clean_%j.err

# 1. Load the same python module you used for generation
# (Check 'module avail python' if 3.9.16 isn't the right one)
module load python/3.9.16s-x3wdtvt

# 2. Activate your environment
source ~/geo_env/bin/activate

# 3. Navigate to your project folder
# Change 'betageometry' to whatever you named your folder
cd ~/betageometry

echo "Starting cleanup of empty JSON files..."

# Run the cleaning script
# -p data : Scan the 'data' folder
# -x      : EXECUTE (actually delete the files)
# -r      : Recursive (optional, safer to include)
python clean_empty_json.py -p data -r -x

echo "Cleanup complete."