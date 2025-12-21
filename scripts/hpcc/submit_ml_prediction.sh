#!/bin/bash
#SBATCH --job-name=ml_predict
#SBATCH --output=logs/ml_prediction_%j.out
#SBATCH --error=logs/ml_prediction_%j.err
#SBATCH --time=02:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --partition=batch

echo "============================================================"
echo "ML Prediction - Script 03"
echo "Job ID: $SLURM_JOB_ID"
echo "Started: $(date)"
echo "============================================================"

# Set project directory (modify as needed)
PROJECT_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
cd ${PROJECT_DIR}
mkdir -p logs

# Activate environment (modify path as needed)
# source activate your_env_name

# Run ML prediction
python scripts/run_03_ml_prediction.py

echo ""
echo "============================================================"
echo "ML Prediction completed: $(date)"
echo "============================================================"

# Check outputs
ls -lh results/tables/ml_*.csv
