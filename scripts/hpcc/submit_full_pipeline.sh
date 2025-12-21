#!/bin/bash
#SBATCH --job-name=target_discovery
#SBATCH --output=logs/full_pipeline_%j.out
#SBATCH --error=logs/full_pipeline_%j.err
#SBATCH --time=24:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=16
#SBATCH --partition=batch

echo "============================================================"
echo "AI-Enabled Therapeutic Target Discovery Pipeline"
echo "Job ID: $SLURM_JOB_ID"
echo "Started: $(date)"
echo "============================================================"

PROJECT_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
cd ${PROJECT_DIR}
mkdir -p logs results/tables figures

# Activate environment (modify as needed)
# source activate target-discovery

# Run pipeline steps
echo "Step 1: Data Preparation"
python scripts/run_01_data_prep.py

echo "Step 2: ML Prediction"
python scripts/run_03_ml_prediction.py

echo "Step 3: Trajectory Analysis"
python scripts/run_05_trajectory.py

echo "Step 4: Target Prioritization"
python scripts/run_06_prioritization.py

echo "Step 5: Drug Repurposing Query"
python scripts/run_08_cmap_query.py

echo ""
echo "============================================================"
echo "Pipeline completed: $(date)"
echo "============================================================"
echo ""
echo "Results saved to:"
ls -lh results/tables/*.csv
echo ""
ls -lh figures/*.png
