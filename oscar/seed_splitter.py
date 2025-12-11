import os
import math

# Configuration
INPUT_FILE = "old_seeds.txt"
OUTPUT_DIR = "seed_batches"
BATCH_SIZE = 36

def main():
    if not os.path.exists(INPUT_FILE):
        print(f"Error: '{INPUT_FILE}' not found. Please upload it first.")
        return

    # 1. Read all valid seeds
    with open(INPUT_FILE, 'r') as f:
        # Filter for non-empty lines
        seeds = [line.strip() for line in f if line.strip()]
    
    total_seeds = len(seeds)
    print(f"Found {total_seeds} seeds in {INPUT_FILE}.")

    # 2. Create the output directory
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # 3. Calculate batches
    # 4577 / 36 = 127.14 -> 128 batches
    num_batches = math.ceil(total_seeds / BATCH_SIZE)
    print(f"Splitting into {num_batches} files (max {BATCH_SIZE} seeds each)...")

    # 4. Write files
    for i in range(num_batches):
        start_idx = i * BATCH_SIZE
        end_idx = start_idx + BATCH_SIZE
        batch_seeds = seeds[start_idx:end_idx]
        
        # Create filename: seed_batches/batch_0.txt, batch_1.txt, etc.
        filename = os.path.join(OUTPUT_DIR, f"batch_{i}.txt")
        
        with open(filename, 'w') as f:
            f.write("\n".join(batch_seeds) + "\n")
            
    print(f"Done. {num_batches} batch files created in '{OUTPUT_DIR}/'.")

if __name__ == "__main__":
    main()