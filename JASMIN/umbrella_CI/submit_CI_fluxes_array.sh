#!/bin/bash

#SBATCH --job-name=CI_fluxes
#SBATCH --account=kscale
#SBATCH --partition=standard
#SBATCH --qos=high
#SBATCH --array=0-5%2
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=23:30:00
#SBATCH --output=/home/users/cornkle/pythonWorkspace/proj_CEH/JASMIN/umbrella_CI/slurm_CI_fluxes_%A_%a.out
#SBATCH --error=/home/users/cornkle/pythonWorkspace/proj_CEH/JASMIN/umbrella_CI/slurm_CI_fluxes_%A_%a.err
#SBATCH --open-mode=append

set -euo pipefail

PROJECT_ROOT="/home/users/cornkle/pythonWorkspace/proj_CEH"
SCRIPT_DIR="${PROJECT_ROOT}/JASMIN/umbrella_CI"

STANDARD_SCRIPT="${SCRIPT_DIR}/umbrella_CI_composites_200km_fluxes.py"
FIXED12_SCRIPT="${SCRIPT_DIR}/umbrella_CI_composites_200km_fluxes_fixed12.py"

source /home/users/cornkle/miniforge3/etc/profile.d/conda.sh
conda activate work

cd "${PROJECT_ROOT}"

export PYTHONUNBUFFERED=1
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

EXPERIMENTS=(standard standard standard fixed12 fixed12 fixed12)
LAGS=(-3 -2 -1 -3 -2 -1)

EXPERIMENT="${EXPERIMENTS[$SLURM_ARRAY_TASK_ID]}"
LAG="${LAGS[$SLURM_ARRAY_TASK_ID]}"

if [[ "${EXPERIMENT}" == "standard" ]]; then
    PYTHON_SCRIPT="${STANDARD_SCRIPT}"
else
    PYTHON_SCRIPT="${FIXED12_SCRIPT}"
fi

echo "================================================================================"
echo "Array job ID : ${SLURM_ARRAY_JOB_ID}"
echo "Array task   : ${SLURM_ARRAY_TASK_ID}"
echo "Node         : $(hostname)"
echo "Experiment   : ${EXPERIMENT}"
echo "Lag          : ${LAG} h"
echo "Script       : ${PYTHON_SCRIPT}"
echo "Started      : $(date --iso-8601=seconds)"
echo "================================================================================"

python -u "${PYTHON_SCRIPT}" --drop-exact-duplicates --lag-hours "${LAG}" --lazy-subsets --checkpoint-every 1 --no-plots

echo "================================================================================"
echo "Experiment   : ${EXPERIMENT}"
echo "Lag          : ${LAG} h"
echo "Completed    : $(date --iso-8601=seconds)"
echo "================================================================================"
