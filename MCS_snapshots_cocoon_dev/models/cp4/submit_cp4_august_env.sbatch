#!/bin/bash
#SBATCH --job-name=cp4_env_aug
#SBATCH --partition=standard
#SBATCH --qos=standard
#SBATCH --account=cp4a
#SBATCH --time=23:45:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --array=0-9%2
#SBATCH --output=logs/env_%A_%a.out
#SBATCH --error=logs/env_%A_%a.err

set -euo pipefail

# IMPORTANT: copy the working --account line from your successful
# submit_cp4_full.sbatch, or supply it at submission: sbatch -A ACCOUNT ...
# This file deliberately does not assume a SLURM account from the storage path.

RUNS=(
 CanESM5-HIST CanESM5-FUT
 HadGEM3-GC31-MM-HIST HadGEM3-GC31-MM-FUT
 MPI-ESM1-2-LR-HIST MPI-ESM1-2-LR-FUT
 MRI-ESM2-0-HIST MRI-ESM2-0-FUT
 TaiESM1-HIST TaiESM1-FUT
)
RUN="${RUNS[$SLURM_ARRAY_TASK_ID]}"

# On JASMIN, BASH_SOURCE points into the Slurm spool, not the project directory.
SCRIPT_DIR="$SLURM_SUBMIT_DIR"
PYTHON="$HOME/miniforge3/envs/work/bin/python"

cd "$SCRIPT_DIR"
echo "Beginning environmental samples $RUN on $(hostname)"
"$PYTHON" "$SCRIPT_DIR/sample_cp4_environments.py" \
    --run "$RUN" --month 8 --all-years \
    --storm-hours 18 --rel-hours -6 0 6 \
    --profile core

echo "Finished $RUN"
