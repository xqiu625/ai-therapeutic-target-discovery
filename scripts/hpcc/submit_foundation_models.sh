#!/bin/bash
#SBATCH --job-name=sepsis_fm
#SBATCH --output=logs/fm_%j.out
#SBATCH --error=logs/fm_%j.err
#SBATCH --time=08:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH --partition=gpu

# Load modules (adjust based on your HPCC)
module purge
module load Conda/3
module load CUDA/11.8

# Activate environment
conda activate sepsis-target

# Create logs directory if not exists
mkdir -p logs

# Check GPU availability
python -c "import torch; print(f'GPU available: {torch.cuda.is_available()}'); print(f'GPU name: {torch.cuda.get_device_name(0) if torch.cuda.is_available() else None}')"

# Run foundation model embeddings notebook
jupyter nbconvert --to notebook --execute \
    --ExecutePreprocessor.timeout=28800 \
    --output notebooks/02_foundation_model_embeddings_executed.ipynb \
    notebooks/02_foundation_model_embeddings.ipynb

echo "Foundation model embeddings completed!"
