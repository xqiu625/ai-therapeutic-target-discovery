#!/bin/bash
#SBATCH --job-name=sepsis_prio
#SBATCH --output=logs/prioritization_%j.out
#SBATCH --error=logs/prioritization_%j.err
#SBATCH --time=02:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --partition=batch

echo "============================================================"
echo "Target Prioritization - Script 06"
echo "Job ID: $SLURM_JOB_ID"
echo "Started: $(date)"
echo "============================================================"

# Set paths
BASE_DIR="/bigdata/godziklab/shared/Xinru/302004"
PROJECT_DIR="${BASE_DIR}/302004_git"

cd ${PROJECT_DIR}
mkdir -p logs

# Activate environment
source activate /bigdata/godziklab/shared/Xinru/pyscenic_env

# Run prioritization
python scripts/run_06_prioritization.py

echo ""
echo "============================================================"
echo "Target prioritization completed: $(date)"
echo "============================================================"

# Check outputs
ls -lh results/tables/prioritized_targets.csv
ls -lh results/tables/top50_targets.csv
