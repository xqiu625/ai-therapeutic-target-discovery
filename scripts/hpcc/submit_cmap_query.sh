#!/bin/bash
#SBATCH --job-name=sepsis_cmap
#SBATCH --output=logs/cmap_query_%j.out
#SBATCH --error=logs/cmap_query_%j.err
#SBATCH --time=01:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2
#SBATCH --partition=batch

echo "============================================================"
echo "CMap/LINCS Drug Query - Script 08"
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

# Run CMap query
python scripts/run_08_cmap_query.py

echo ""
echo "============================================================"
echo "CMap query completed: $(date)"
echo "============================================================"

# Check outputs
ls -lh results/tables/cmap_*.csv
ls -lh results/tables/sigcom_*.csv
