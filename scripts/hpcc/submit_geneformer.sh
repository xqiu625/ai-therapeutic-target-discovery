#!/bin/bash
#SBATCH --job-name=sepsis_geneformer
#SBATCH --output=logs/geneformer_%j.out
#SBATCH --error=logs/geneformer_%j.err
#SBATCH --time=04:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH --partition=gpu

# ============================================================
# Geneformer Embedding Generation
# Environment: Singularity container
# ============================================================

echo "============================================================"
echo "Geneformer Embedding Generation"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo "Started: $(date)"
echo "============================================================"

# Set paths
BASE_DIR="/bigdata/godziklab/shared/Xinru/302004"
PROJECT_DIR="${BASE_DIR}/302004_git"
INPUT_FILE="${PROJECT_DIR}/data/processed/adata_processed.h5ad"
OUTPUT_DIR="${PROJECT_DIR}/data/embeddings"

# Geneformer container path
GENEFORMER_SIF="/bigdata/godziklab/shared/Xinru/302006/07_CONTAINERS/geneformer.sif"

# Create directories
mkdir -p ${PROJECT_DIR}/logs
mkdir -p ${OUTPUT_DIR}

# Load modules
module purge
module load Singularity

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

# Check if container exists
if [ ! -f "$GENEFORMER_SIF" ]; then
    echo "ERROR: Geneformer container not found: $GENEFORMER_SIF"
    exit 1
fi

# Run Geneformer via Singularity
cd ${PROJECT_DIR}

echo ""
echo "Running Geneformer embedding generation..."
singularity exec --nv \
    --bind ${BASE_DIR}:${BASE_DIR} \
    ${GENEFORMER_SIF} \
    python scripts/run_geneformer.py \
        --input ${INPUT_FILE} \
        --output_dir ${OUTPUT_DIR} \
        --model_dir ctheodoris/Geneformer \
        --batch_size 32

echo ""
echo "============================================================"
echo "Geneformer completed: $(date)"
echo "============================================================"
