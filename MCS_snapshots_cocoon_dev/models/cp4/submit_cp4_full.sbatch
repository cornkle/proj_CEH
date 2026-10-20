#!/bin/bash
#SBATCH --job-name=cp4_cocoon
#SBATCH --partition=standard
#SBATCH --qos=standard
#SBATCH --account=cp4a
#SBATCH --time=23:45:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --array=0-9%3
#SBATCH --output=logs/cp4_%A_%a.out
#SBATCH --error=logs/cp4_%A_%a.err

set -euo pipefail

# Array indices are deliberately a fixed model/experiment mapping.
RUNS=(
    CanESM5-HIST
    CanESM5-FUT
    HadGEM3-GC31-MM-HIST
    HadGEM3-GC31-MM-FUT
    MPI-ESM1-2-LR-HIST
    MPI-ESM1-2-LR-FUT
    MRI-ESM2-0-HIST
    MRI-ESM2-0-FUT
    TaiESM1-HIST
    TaiESM1-FUT
)

RUN="${RUNS[$SLURM_ARRAY_TASK_ID]}"
SCRIPT_DIR="$SLURM_SUBMIT_DIR"
PYTHON="$HOME/miniforge3/envs/work/bin/python"

if [[ ! -x "$PYTHON" ]]; then
    echo "Python not found: $PYTHON" >&2
    echo "Edit PYTHON in this file to match your environment." >&2
    exit 1
fi

echo "Started $(date -u) UTC | run=$RUN | job=${SLURM_JOB_ID:-unknown}"
echo "Python: $PYTHON"
echo "Script: $SCRIPT_DIR/build_cp4_base_tables.py"

# First runs August in every available year, then all other months.
# Already completed monthly Parquets are automatically skipped.
"$PYTHON" "$SCRIPT_DIR/build_cp4_base_tables.py" --run "$RUN" --all

echo "Finished $(date -u) UTC | run=$RUN"
