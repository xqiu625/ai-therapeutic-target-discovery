#!/usr/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=10 
#SBATCH --mem=256gbG
#SBATCH --job-name="inter_DC"
#SBATCH -p intel,batch,short 

module unload R/4.0.1
module load R/4.1.0_gcc-8.3.0
module unload Rcpp
module load hdf5/1.12.0_gcc-8.3.0
module load workspace/scratch
export TMPDIR=$SCRATCH

Rscript --vanilla cell-cell_interaction_score_by_Sample_DC.R