"""Split a seed range into per-task batch files for a SLURM array job.

    python split_seeds.py --range 0-9999 --tasks 128
    python split_seeds.py --from-file seeds.txt --batch-size 36

Writes `seed_batches/batch_0.txt` ... and prints the `--array` range to use.
"""

import argparse
import math
import os
import sys


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    src = parser.add_mutually_exclusive_group(required=True)
    src.add_argument("--range", dest="seed_range", help="inclusive range, e.g. 0-9999")
    src.add_argument("--from-file", help="file with one seed per line")
    size = parser.add_mutually_exclusive_group()
    size.add_argument("--tasks", type=int, help="number of array tasks to split across")
    size.add_argument("--batch-size", type=int, help="seeds per task")
    parser.add_argument("--out-dir", default="seed_batches")
    args = parser.parse_args()

    if args.seed_range:
        lo, _, hi = args.seed_range.partition("-")
        seeds = [str(s) for s in range(int(lo), int(hi) + 1)]
    else:
        if not os.path.exists(args.from_file):
            sys.exit(f"Not found: {args.from_file}")
        with open(args.from_file, "r", encoding="utf-8") as f:
            seeds = [line.strip() for line in f if line.strip()]

    if not seeds:
        sys.exit("No seeds to split.")

    if args.tasks:
        batch_size = math.ceil(len(seeds) / args.tasks)
    elif args.batch_size:
        batch_size = args.batch_size
    else:
        batch_size = math.ceil(len(seeds) / 100)

    os.makedirs(args.out_dir, exist_ok=True)
    num_batches = math.ceil(len(seeds) / batch_size)
    for i in range(num_batches):
        chunk = seeds[i * batch_size:(i + 1) * batch_size]
        with open(os.path.join(args.out_dir, f"batch_{i}.txt"), "w", encoding="utf-8") as f:
            f.write("\n".join(chunk) + "\n")

    print(f"{len(seeds)} seeds -> {num_batches} files of up to {batch_size} in {args.out_dir}/")
    print(f"Use:  #SBATCH --array=0-{num_batches - 1}")


if __name__ == "__main__":
    main()
