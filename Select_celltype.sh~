#!/usr/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=10 
#SBATCH --mem=512gbG
#SBATCH --job-name="SeuV4D7"
#SBATCH --time=4-00:00:00
#SBATCH -p highmem

module unload R/4.0.1
module load R/4.1.0_gcc-8.3.0
module unload Rcpp
module load hdf5/1.12.0_gcc-8.3.0
module load workspace/scratch
export TMPDIR=$SCRATCH

Rscript --vanilla SeuratV4_celltypes_dat7.R