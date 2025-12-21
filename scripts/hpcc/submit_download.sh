#!/bin/bash
#SBATCH --job-name=sepsis_download
#SBATCH --output=logs/download_%j.out
#SBATCH --error=logs/download_%j.err
#SBATCH --time=04:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2
#SBATCH --partition=standard

# Load modules
module purge
module load Conda/3

# Activate environment
conda activate sepsis-target

# Create directories
mkdir -p logs data/raw data/external

# Download GEO data (non-interactive mode)
echo "Downloading GSE167363 data..."
cd data/raw

# Download the RAW tar file
wget -c "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE167nnn/GSE167363/suppl/GSE167363_RAW.tar"

# Extract
tar -xvf GSE167363_RAW.tar

cd ../..

# Download pySCENIC databases
echo "Downloading pySCENIC databases..."
cd data/external

# TF list
wget -c "https://raw.githubusercontent.com/aertslab/pySCENIC/master/resources/hs_hgnc_tfs.txt"

# cisTarget databases (these are large)
wget -c "https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc_v10_clust/gene_based/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"
wget -c "https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc_v10_clust/gene_based/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"
wget -c "https://resources.aertslab.org/cistarget/motif2tf/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl"

cd ../..

echo "Download complete!"
ls -la data/raw/
ls -la data/external/
