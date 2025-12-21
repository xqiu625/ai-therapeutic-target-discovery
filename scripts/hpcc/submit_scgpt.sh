#!/bin/bash
#SBATCH --job-name=sepsis_scgpt
#SBATCH --output=logs/scgpt_%j.out
#SBATCH --error=logs/scgpt_%j.err
#SBATCH --time=04:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH --partition=gpu

# ============================================================
# scGPT Embedding Generation
# Environment: scgpt (conda)
# ============================================================

echo "============================================================"
echo "scGPT Embedding Generation"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo "Started: $(date)"
echo "============================================================"

# Set paths
BASE_DIR="/bigdata/godziklab/shared/Xinru/302004"
PROJECT_DIR="${BASE_DIR}/302004_git"
INPUT_FILE="${PROJECT_DIR}/data/processed/adata_processed.h5ad"
OUTPUT_DIR="${PROJECT_DIR}/data/embeddings"

# Create directories
mkdir -p ${PROJECT_DIR}/logs
mkdir -p ${OUTPUT_DIR}

# Load modules (adjust as needed for your HPC)
module purge
module load Conda/3
module load CUDA/11.8

# Activate scGPT environment
conda activate scgpt

# Check GPU
echo ""
echo "GPU Check:"
nvidia-smi

# Check if input file exists
if [ ! -f "$INPUT_FILE" ]; then
    echo "ERROR: Input file not found: $INPUT_FILE"
    echo "Please run notebook 01 first to create the processed data."
    exit 1
fi

# Run scGPT
cd ${PROJECT_DIR}

echo ""
echo "Running scGPT embedding generation..."
python scripts/run_scgpt.py \
    --input ${INPUT_FILE} \
    --output_dir ${OUTPUT_DIR} \
    --model_dir scGPT_human \
    --batch_size 64

echo ""
echo "============================================================"
echo "scGPT completed: $(date)"
echo "============================================================"
