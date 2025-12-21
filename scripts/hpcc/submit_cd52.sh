#!/bin/bash
#SBATCH --job-name=sepsis_cd52
#SBATCH --output=logs/cd52_%j.out
#SBATCH --error=logs/cd52_%j.err
#SBATCH --time=02:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --partition=batch

echo "============================================================"
echo "CD52 Deep Dive - Script 07"
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

# Run CD52 deep dive
python scripts/run_07_cd52_deep_dive.py

echo ""
echo "============================================================"
echo "CD52 deep dive completed: $(date)"
echo "============================================================"

# Check outputs
ls -lh figures/fig7_cd52_summary.png
ls -lh results/tables/cd52_*.csv
