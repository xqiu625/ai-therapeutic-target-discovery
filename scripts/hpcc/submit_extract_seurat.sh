#!/bin/bash
#SBATCH --job-name=extract_seurat
#SBATCH --output=logs/extract_seurat_%j.out
#SBATCH --error=logs/extract_seurat_%j.err
#SBATCH --time=04:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=4
#SBATCH --partition=batch

# ============================================================
# Extract Seurat Data to CSV/MTX for Python
# ============================================================

echo "============================================================"
echo "Seurat Data Extraction"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo "Started: $(date)"
echo "============================================================"

# Set paths
BASE_DIR="/bigdata/godziklab/shared/Xinru/302004"
PROJECT_DIR="${BASE_DIR}/302004_git"

# Create directories
mkdir -p ${PROJECT_DIR}/logs
mkdir -p ${PROJECT_DIR}/data/processed

# Load R module (adjust based on your HPC)
module purge
module load R

# Check R version
echo ""
echo "R version:"
R --version | head -1

# Run extraction script
cd ${PROJECT_DIR}

echo ""
echo "Running Seurat extraction..."
echo ""

Rscript scripts/extract_seurat_data.R

echo ""
echo "============================================================"
echo "Extraction completed: $(date)"
echo "============================================================"

# List output files
echo ""
echo "Output files:"
ls -lh ${PROJECT_DIR}/data/processed/
