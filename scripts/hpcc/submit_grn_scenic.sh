#!/bin/bash
#SBATCH --job-name=sepsis_scenic
#SBATCH --output=logs/grn_scenic_%j.out
#SBATCH --error=logs/grn_scenic_%j.err
#SBATCH --time=24:00:00
#SBATCH --mem=256G
#SBATCH --cpus-per-task=16
#SBATCH --partition=batch

echo "============================================================"
echo "SCENIC GRN Analysis - Script 04"
echo "Job ID: $SLURM_JOB_ID"
echo "Started: $(date)"
echo "============================================================"

# Set paths
BASE_DIR="/bigdata/godziklab/shared/Xinru/302004"
PROJECT_DIR="${BASE_DIR}/302004_git"

cd ${PROJECT_DIR}
mkdir -p logs
mkdir -p data/grn
mkdir -p data/external

# Activate environment
source activate /bigdata/godziklab/shared/Xinru/pyscenic_env

# Check for ranking databases (user needs to download manually if missing)
echo ""
echo "Checking SCENIC resources..."
DB_COUNT=$(ls ${PROJECT_DIR}/data/external/*.feather 2>/dev/null | wc -l)
if [ "$DB_COUNT" -eq 0 ]; then
    echo ""
    echo "WARNING: No ranking databases found!"
    echo "SCENIC cisTarget step will be skipped."
    echo "To enable, download from: https://resources.aertslab.org/cistarget/databases/"
    echo ""
fi

# Run SCENIC analysis
python scripts/run_04_grn_scenic.py

echo ""
echo "============================================================"
echo "SCENIC analysis completed: $(date)"
echo "============================================================"

# Check outputs
ls -lh data/grn/
ls -lh results/tables/
