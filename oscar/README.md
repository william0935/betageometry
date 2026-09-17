# Running BetaGeometry on OSCAR

Generating training data and fine-tuning Gemma on Brown's OSCAR cluster.

## 1. Prerequisites

An OSCAR account — [Brown CCV docs](https://docs.ccv.brown.edu/oscar).

## 2. Upload

Use [Open OnDemand](https://ood.ccv.brown.edu) → **Files** → **Home Directory**.

Create `~/betageometry` and upload the repository into it, **keeping the directory
structure** — `oscar/*.sh` reference `oscar/` paths and the Python modules import each
other by name from the repository root:

```text
~/betageometry/
├── generate_batch.py
├── consolidate_data.py
├── split_seeds.py
├── finetune.py
├── data_generation.py
├── constructions.py  problem.py  dd_ar.py  ar.py  relations.py  rabbits.py  gemma.py
└── oscar/
    ├── submit_generate.sh
    ├── submit_consolidate.sh
    └── submit_finetune.sh
```

## 3. Environment (first run only)

```bash
module avail python           # confirm the version below still exists
module load python/3.11.0s-ixrhc3q

# The virtual environment lives outside the repo so re-uploading the code cannot
# clobber it, and so every job can share one installation.
python3 -m venv ~/geo_env
source ~/geo_env/bin/activate

cd ~/betageometry
python -m pip install -e .            # numpy + matplotlib, enough for data generation
python -m pip install -e '.[llm]'     # adds torch/transformers/peft, for fine-tuning only
```

Each `submit_*.sh` activates `~/geo_env` itself, so the batch jobs do not depend on
your shell. Reactivate it by hand whenever you run something interactively:

```bash
source ~/geo_env/bin/activate
```

## 4. Generate data

Split the seed range into one file per array task. The command prints the `--array`
range to use:

```bash
cd ~/betageometry
python split_seeds.py --range 0-99999 --tasks 128
# -> 100000 seeds -> 128 files of up to 782 in seed_batches/
# -> Use:  #SBATCH --array=0-127
```

Set `--array` in `oscar/submit_generate.sh` to match, then submit:

```bash
sbatch oscar/submit_generate.sh     # note the job ID it prints
```

Each task writes `data/training_data_<seed>.json` for every seed that produced
examples, and `data/productive_batch_<n>.txt` listing those seeds. Seeds that produce
nothing write no file, so there is nothing empty to clean up afterwards.

## 5. Consolidate

Chain it after the array job so it waits:

```bash
sbatch --dependency=afterok:<ARRAY_JOB_ID> oscar/submit_consolidate.sh
```

This merges every per-seed file into `data/training_data.json`, dropping duplicate
(input, output) pairs. Add `--clean` to `consolidate_data.py` in the script if you want
the per-seed files deleted after a successful merge.

Check the result before training on it:

```bash
python -c "import json; d=json.load(open('data/training_data.json')); print(len(d), 'examples'); print(d[0])"
```

## 6. Fine-tune

Gemma is gated: accept the licence on HuggingFace, then authenticate once on the
cluster.

```bash
huggingface-cli login
sbatch oscar/submit_finetune.sh
```

The LoRA adapter lands in `gemma-finetuned-geometry/`. Point the solver at it:

```bash
python solve.py problem1 --gemma gemma-finetuned-geometry
```

## 7. Download

**Files** → `~/betageometry/data/` → select `training_data.json` → **Download**.

---

### Troubleshooting

**`Batch script contains DOS line breaks (\r\n)`** — uploaded from Windows:

```bash
dos2unix oscar/*.sh
```

**Jobs die at import with a display error** — a compute node has no X display. The
entry points already call `matplotlib.use("Agg")` before importing pyplot; if you add a
new script, do the same, before any other project import.

**`CUDA out of memory` during fine-tuning** — lower `--batch-size` and raise
`--grad-accumulation` to keep the effective batch size constant.
