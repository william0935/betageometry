"""Generate training data for a batch of seeds (the Oscar array-job entry point).

Each SLURM array task is handed a file of seeds and writes one JSON file per seed into
`data/`, plus a summary listing the seeds that produced anything. `consolidate_data.py`
merges them afterwards.

    python generate_batch.py seed_batches/batch_0.txt
    python generate_batch.py --seeds 0-999            # contiguous range, no batch file
"""

import argparse
import os
import sys
import traceback
from typing import List

import matplotlib
matplotlib.use("Agg")  # must precede any pyplot import; compute nodes have no display

from data_generation import generate_for_seed, write_examples


def read_seed_file(path: str) -> List[int]:
    seeds = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if line.isdigit():
                seeds.append(int(line))
    return seeds


def parse_range(spec: str) -> List[int]:
    if "-" in spec:
        lo, hi = spec.split("-", 1)
        return list(range(int(lo), int(hi) + 1))
    return [int(spec)]


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("seed_file", nargs="?", help="file with one seed per line")
    parser.add_argument("--seeds", help="inclusive range like 0-999, instead of a file")
    parser.add_argument("--rounds", type=int, default=4,
                        help="auxiliary points added per configuration")
    parser.add_argument("--data-dir", default="data")
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args()

    if args.seeds:
        seeds = parse_range(args.seeds)
        batch_name = f"range_{args.seeds}"
    elif args.seed_file:
        if not os.path.exists(args.seed_file):
            sys.exit(f"Seed file not found: {args.seed_file}")
        seeds = read_seed_file(args.seed_file)
        batch_name = os.path.splitext(os.path.basename(args.seed_file))[0]
    else:
        sys.exit("Provide a seed file or --seeds RANGE")

    os.makedirs(args.data_dir, exist_ok=True)
    print(f"Processing {len(seeds)} seeds.")

    productive = []
    total = 0
    for seed in seeds:
        try:
            examples = generate_for_seed(seed, rounds=args.rounds, verbose=False)
        except Exception:
            # One bad configuration must not take down a job holding hundreds of seeds.
            print(f"seed {seed}: FAILED", file=sys.stderr)
            traceback.print_exc()
            continue

        if not examples:
            if not args.quiet:
                print(f"seed {seed}: 0 examples")
            continue

        # Only seeds that produced something get a file, so there is nothing empty for
        # the consolidation step to clean up afterwards.
        write_examples(examples, os.path.join(args.data_dir, f"training_data_{seed}.json"))
        productive.append(seed)
        total += len(examples)
        if not args.quiet:
            print(f"seed {seed}: {len(examples)} examples")

    summary = os.path.join(args.data_dir, f"productive_{batch_name}.txt")
    with open(summary, "w", encoding="utf-8") as f:
        f.write("\n".join(str(s) for s in productive) + ("\n" if productive else ""))

    print(f"Done: {total} examples from {len(productive)}/{len(seeds)} seeds.")
    print(f"Productive seeds listed in {summary}")


if __name__ == "__main__":
    main()
