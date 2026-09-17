<<<<<<< HEAD
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

python3 -m venv ~/geo_env
source ~/geo_env/bin/activate

cd ~/betageometry
pip install -e .              # numpy + matplotlib, enough for data generation
pip install -e '.[llm]'       # adds torch/transformers/peft — only needed for fine-tuning
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
=======
# Synthetic Geometry Data Generation on OSCAR

This repository contains the scripts necessary to generate synthetic training data for geometry problems using Brown University's OSCAR high-performance computing cluster.

## 1\. Prerequisites

Before running these scripts, you must have an active account on the OSCAR cluster.

  * **Request an account/Documentation:** [Brown CCV Documentation](https://docs.ccv.brown.edu/oscar)

## 2\. Access & File Upload

We recommend using **Open OnDemand (OOD)** for easy file management and terminal access.

1.  Log in to [https://ood.ccv.brown.edu](https://ood.ccv.brown.edu).
2.  In the top menu, navigate to **Files** \> **Home Directory**.
3.  Create a new folder named `betageometry`.
4.  **Upload the files:**
      * Navigate inside `betageometry`.
      * Upload all Python scripts (`.py`) and Shell scripts (`.sh`) to this repository.
      * **Important:** Upload the files directly into `betageometry`. Do **not** create a subfolder named `oscar/` or `src/`.

**Correct File Structure:**

```text
~/betageometry/
├── oscar_data_generation.py
├── training_file_generation.py
├── submit_array.sh
├── submit_consolidate.sh
├── Constructions.py
├── Problem.py
└── ... (other helper modules)
```

## 3\. Environment Setup (First Run Only)

Before submitting jobs, you must set up the Python environment.

1.  In Open OnDemand, go to **Clusters** \> **Shell Access**.
2.  Run the following commands to create the virtual environment:

<!-- end list -->

```bash
# Load Python (check 'module avail python' for specific versions)
module load python/3.9.16s-x3wdtvt

# Create virtual environment
python3 -m venv ~/geo_env

# Activate and install dependencies
source ~/geo_env/bin/activate
pip install numpy matplotlib
```

## 4\. Running the Generation (Parallel)

The data generation is split into parallel "array" jobs to speed up the process.

1.  Navigate to your folder:
    ```bash
    cd ~/betageometry
    ```
2.  **Submit the Job Array:**
    This will launch multiple workers (e.g., 100 tasks) to process seeds in parallel.
    ```bash
    sbatch submit_array.sh
    ```
    *Output:* `Submitted batch job 12345678` (Note this Job ID).

## 5\. Consolidating the Data

Once the generation jobs finish, you need to merge the thousands of small JSON files into one training file.

**Option A: The Automatic Way (Recommended)**
Submit this immediately after the array job, using the Job ID you just got. It will wait in the queue until the generation is 100% complete.

```bash
# Replace 12345678 with your actual Array Job ID
sbatch --dependency=afterok:12345678 submit_consolidate.sh
```

**Option B: The Manual Way**
Wait for the array jobs to finish (check using `squeue -u $USER`), then download the whole `data` folder to run `training_file_generation.py`.

## 6\. Downloading Results

Once the consolidation job is complete:

1.  Go back to the **Open OnDemand Files** dashboard.
2.  Navigate to `~/betageometry/data/`.
3.  Locate `training_data.json` (this is the final dataset).
4.  Select the `data` folder or the specific JSON file and click **Download**. Downloading the whole `data` folder is recommended.

-----

### Troubleshooting

**Error: `Batch script contains DOS line breaks (\r\n)`**
If you uploaded files from Windows, they might have the wrong line endings. Fix them on the cluster by running:

```bash
dos2unix submit_array.sh submit_consolidate.sh
```
>>>>>>> 6b5bb16886069cf2d557723cbaee5eb3172e2f42
