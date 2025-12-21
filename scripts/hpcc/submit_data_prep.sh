#!/bin/bash
#SBATCH --job-name=sepsis_data_prep
#SBATCH --output=logs/data_prep_%j.out
#SBATCH --error=logs/data_prep_%j.err
#SBATCH --time=04:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=4
#SBATCH --partition=batch

echo "============================================================"
echo "Data Preparation - Notebook 01"
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

# Run data preparation script
python scripts/run_01_data_prep.py

echo ""
echo "============================================================"
echo "Data preparation completed: $(date)"
echo "============================================================"

# Check output
ls -lh data/processed/
