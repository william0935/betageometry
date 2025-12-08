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