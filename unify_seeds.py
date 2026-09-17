from pathlib import Path
import re
import argparse
from collections import OrderedDict

TOKEN_SPLIT_RE = re.compile(r"[,\s;]+")


def find_seed_files(root: Path, exclude: Path = None):
    root = root.resolve()
    for p in root.rglob("*"):
        if not p.is_file():
            continue
        name = p.name.lower()
        if "nontrivial_seeds" in name:
            if exclude and p.resolve() == exclude.resolve():
                continue
            yield p


def parse_seeds_from_file(path: Path):
    seeds = []
    text = path.read_text(encoding="utf-8", errors="ignore")
    for raw_line in text.splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith("#") or line.startswith("//"):
            continue
        # split tokens by common delimiters, keep non-empty tokens
        for tok in TOKEN_SPLIT_RE.split(line):
            tok = tok.strip()
            if not tok:
                continue
            seeds.append(tok)
    return seeds


def collect_seeds(files):
    seen = OrderedDict()
    for f in files:
        for s in parse_seeds_from_file(f):
            if s not in seen:
                seen[s] = None
    return list(seen.keys())


def main():
    parser = argparse.ArgumentParser(description="Collect seed values into one seeds.txt")
    parser.add_argument(
        "--out",
        "-o",
        type=Path,
        default=None,
        help="Output file path (defaults to project_root/seeds.txt)",
    )
    args = parser.parse_args()

    script_dir = Path(__file__).parent.resolve()  # data/
    project_root = script_dir.parent.resolve()
    out_path = args.out.resolve() if args.out else (project_root / "seeds.txt")
    # Support running this script either from the project root or from the data/ subdirectory.
    if script_dir.name.lower() == "data":
        project_root = script_dir.parent.resolve()
    else:
        project_root = script_dir.resolve()

    # Recompute out_path when not explicitly provided
    out_path = args.out.resolve() if args.out else (project_root / "seeds.txt")
    seed_files = list(find_seed_files(project_root, exclude=out_path))
    if not seed_files:
        print("No seed files found.")
        return

    seeds = collect_seeds(seed_files)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text("\n".join(seeds) + ("\n" if seeds else ""), encoding="utf-8")

    print(f"Wrote {len(seeds)} unique seeds to: {out_path}")
    print("Found seed files:")
    for f in seed_files:
        print(" -", f.relative_to(project_root))


if __name__ == "__main__":
    main()