#!/bin/bash
#SBATCH --job-name=sepsis_traj
#SBATCH --output=logs/trajectory_%j.out
#SBATCH --error=logs/trajectory_%j.err
#SBATCH --time=04:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --partition=batch

echo "============================================================"
echo "Trajectory Analysis - Script 05"
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

# Run trajectory analysis
python scripts/run_05_trajectory.py

echo ""
echo "============================================================"
echo "Trajectory analysis completed: $(date)"
echo "============================================================"

# Check outputs
ls -lh data/processed/adata_trajectory.h5ad
ls -lh results/tables/trajectory_*.csv
