#!/bin/bash

#SBATCH --mail-user=laura.l.colbran@vanderbilt.edu
#SBATCH --mail-type=FAIL
#SBATCH --time=16:00:00
#SBATCH --nodes=1
#SBATCH --mem=22G
#SBATCH --array=1-22%5
#SBATCH --output=%x_%A_%a.out

#################
# Job Execution #
#################
source /etc/profile.d/lmod.sh

chr=$SLURM_ARRAY_TASK_ID
genes=$gene_annot$SLURM_ARRAY_TASK_ID'.txt'
snps=$snp_annot$SLURM_ARRAY_TASK_ID'.RDS'

module load GCC/5.4.0-2.26  OpenMPI/1.10.3 R/3.3.3
Rscript ./create_model.R $study $expr_RDS $geno $genes $snps \
    $n_k_folds $alpha $out_dir $chr $snpset $window
