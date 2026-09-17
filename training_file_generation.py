import os
import json
import glob

DATA_DIR = os.path.join(os.path.dirname(__file__), "data")
TRAINING_FILE = os.path.join(DATA_DIR, "training_data.json")


def load_existing_training():
    if not os.path.exists(TRAINING_FILE):
        return []
    try:
        with open(TRAINING_FILE, "r", encoding="utf-8") as f:
            data = json.load(f)
    except Exception:
        # Corrupt or unreadable file; treat as empty
        return []
    if isinstance(data, list):
        return data
    if isinstance(data, dict):
        return [data]
    return []


def read_seed_file(path):
    seeds = set()
    try:
        # try JSON first (in case the seed file is a JSON array)
        with open(path, "r", encoding="utf-8") as f:
            text = f.read().strip()
            if not text:
                return seeds
            try:
                parsed = json.loads(text)
                if isinstance(parsed, (list, tuple)):
                    for s in parsed:
                        seeds.add(str(s).strip())
                    return seeds
            except Exception:
                # not JSON array; fall back to line-based parsing
                pass
        # line-based parsing
        with open(path, "r", encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                # allow comma-separated lines
                parts = [p.strip() for p in line.replace(",", " ").split()]
                for p in parts:
                    if p:
                        seeds.add(p)
    except Exception:
        pass
    return seeds


def find_jsons_for_seed(seed):
    # find any json file in DATA_DIR whose filename contains the seed
    pattern = os.path.join(DATA_DIR, f"*{seed}*.json")
    return glob.glob(pattern)


def load_json_file(path):
    try:
        with open(path, "r", encoding="utf-8") as f:
            return json.load(f)
    except Exception:
        return None


def main():
    os.makedirs(DATA_DIR, exist_ok=True)

    # existing = load_existing_training()
    # alternatively, start fresh
    existing = []

    # gather all seeds from the file seeds.txt
    seed_file_path = os.path.join(DATA_DIR, "seeds.txt")
    all_seeds = read_seed_file(seed_file_path)
    print(f"Found {len(all_seeds)} unique seeds from '{seed_file_path}'.")

    # for each seed find matching json files and collect their data
    new_entries = []
    for seed in sorted(all_seeds):
        json_paths = find_jsons_for_seed(seed)
        if not json_paths:
            # no matching json for this seed; skip quietly
            continue
        for jp in json_paths:
            content = load_json_file(jp)
            if content is None:
                continue
            if isinstance(content, list):
                new_entries.extend(content)
            elif isinstance(content, dict):
                new_entries.append(content)
            else:
                # ignore unknown types
                continue

    # merge with existing data (append)
    combined = list(existing) + new_entries

    # write back
    try:
        with open(TRAINING_FILE, "w", encoding="utf-8") as f:
            json.dump(combined, f, indent=2, ensure_ascii=False)
    except Exception as e:
        print(f"Failed to write training file: {e}")
        return

    print(f"{len(combined)} entries in training data.")


if __name__ == "__main__":
    main()