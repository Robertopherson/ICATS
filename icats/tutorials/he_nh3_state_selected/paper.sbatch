#!/bin/bash
#SBATCH --job-name=nh3_icats_tutorial
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=4G
#SBATCH --time=06:00:00
#SBATCH --array=0-255%4
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export PYTHONDONTWRITEBYTECODE=1
"${ICATS_PYTHON:-python3}" -u run.py --paper-task "$SLURM_ARRAY_TASK_ID" \
    --workers "$SLURM_CPUS_PER_TASK" --output paper
