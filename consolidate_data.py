"""Merge the per-seed JSON files from a generation run into one training file.

Replaces BetaGeometry's `training_file_generation.py` + `clean_empty_json.py` pair. That
version globbed `*{seed}*.json` per seed, which matches by substring -- seed 7 also picks
up `training_data_7.json`, `training_data_17.json` and `training_data_700.json` -- so
files were pulled in repeatedly and the merged set was inflated with duplicates. This
walks the directory once instead.

    python consolidate_data.py                    # data/*.json -> data/training_data.json
    python consolidate_data.py --clean            # also delete the per-seed files
"""

import argparse
import glob
import json
import os
import sys
from typing import Dict, List

DEFAULT_DATA_DIR = "data"
OUTPUT_NAME = "training_data.json"


def load_json(path: str):
    try:
        with open(path, "r", encoding="utf-8") as f:
            return json.load(f)
    except Exception as exc:
        print(f"  skipping {path}: {exc}", file=sys.stderr)
        return None


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--data-dir", default=DEFAULT_DATA_DIR)
    parser.add_argument("--out", default=None)
    parser.add_argument("--clean", action="store_true",
                        help="delete the per-seed files after a successful merge")
    parser.add_argument("--keep-duplicates", action="store_true",
                        help="keep identical (input, output) pairs instead of collapsing them")
    args = parser.parse_args()

    out_path = args.out or os.path.join(args.data_dir, OUTPUT_NAME)
    if not os.path.isdir(args.data_dir):
        sys.exit(f"Not a directory: {args.data_dir}")

    sources = sorted(
        p for p in glob.glob(os.path.join(args.data_dir, "*.json"))
        if os.path.abspath(p) != os.path.abspath(out_path)
    )
    print(f"Merging {len(sources)} files from {args.data_dir}/ ...")

    merged: List[Dict[str, str]] = []
    seen = set()
    empty_files = []
    for path in sources:
        content = load_json(path)
        if content is None:
            continue
        if isinstance(content, dict):
            content = [content]
        if not isinstance(content, list) or not content:
            empty_files.append(path)
            continue
        for item in content:
            if not isinstance(item, dict) or "input_text" not in item:
                continue
            if args.keep_duplicates:
                merged.append(item)
                continue
            key = (item["input_text"], item.get("output_text", ""))
            if key in seen:
                continue
            seen.add(key)
            merged.append(item)

    os.makedirs(args.data_dir, exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(merged, f, indent=2, ensure_ascii=False)

    print(f"{len(merged)} examples -> {out_path}")
    if empty_files:
        print(f"{len(empty_files)} source files were empty.")

    if args.clean:
        removed = 0
        for path in sources:
            try:
                os.remove(path)
                removed += 1
            except OSError as exc:
                print(f"  could not remove {path}: {exc}", file=sys.stderr)
        print(f"Removed {removed} per-seed files.")


if __name__ == "__main__":
    main()
