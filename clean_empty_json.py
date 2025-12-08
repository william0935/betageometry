"""
clean_empty_json.py

Remove JSON files that contain an empty object ({}) in a directory.

Usage:
    python clean_empty_json.py         # dry-run, shows files that would be removed
    python clean_empty_json.py -x      # actually deletes matching files
    python clean_empty_json.py -p data -r -x
"""

from pathlib import Path
import argparse
import json
import sys

def is_empty_json_file(path: Path) -> bool:
    try:
        text = path.read_text(encoding="utf-8").strip()
    except Exception:
        return False
    # Quick check for the common case
    if text == "[]":
        return True
        # Fall back to JSON parsing for formatted/whitespace variants
    try:
        obj = json.loads(text)
    except Exception:
        return False
    return isinstance(obj, dict) and len(obj) == 0

def find_json_files(directory: Path, recursive: bool):
    if recursive:
        yield from directory.rglob("*.json")
    else:
        yield from directory.glob("*.json")

def main():
    p = argparse.ArgumentParser(description="Delete JSON files that are empty objects ({})")
    p.add_argument("-p", "--path", default="data", help="Directory to scan (default: data)")
    p.add_argument("-r", "--recursive", action="store_true", help="Recurse into subdirectories")
    p.add_argument("-x", "--execute", action="store_true", help="Actually delete files (default is dry-run)")
    args = p.parse_args()

    base = Path(args.path)
    if not base.exists() or not base.is_dir():
        print(f"Error: path not found or not a directory: {base}", file=sys.stderr)
        sys.exit(2)

    removed = 0
    checked = 0
    for fp in find_json_files(base, args.recursive):
        checked += 1
        if is_empty_json_file(fp):
            if args.execute:
                try:
                    fp.unlink()
                    print(f"Removed: {fp}")
                    removed += 1
                except Exception as e:
                    print(f"Failed to remove {fp}: {e}", file=sys.stderr)
            else:
                print(f"Would remove: {fp}")

    summary = f"Scanned {checked} .json files. " + (f"Removed {removed} empty json files." if args.execute else "Dry-run complete.")
    print(summary)

if __name__ == "__main__":
        main()