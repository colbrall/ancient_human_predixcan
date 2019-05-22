#!/bin/bash
###############################
# Resource Manager Directives #
###############################

### Run using command -v tissue=tissue_name, chrom=chromosome_number, alpha=alpha, window=window

### Select the shell you would like to script to execute within.
##PBS -S /bin/bash
#SBATCH --mail-user=laura.l.colbran@vanderbilt.edu
#SBATCH --mail-type=ALL
### Inform the scheduler of the expected runtime, where
### walltime=HH:MM:SS.
#SBATCH --time=16:00:00
### Inform the scheduler of the number of CPU cores for your job.
#SBATCH --nodes=1
### Inform the scheduler of the amount of memory you expect
### to use. Use units of ‘b’, ‘kb’, ‘mb’, or ‘gb’.
#SBATCH --mem=22G
### Set the destination for your program’s output: stdout and stderr.
##PBS -e ./${PBS_JOBNAME}.e${PBS_JOBID}
##SBATCH --output ./${SLURM_JOB_NAME}.${SLURM_JOB_ID}.out
 
#################
# Job Execution #
#################

module load GCC/5.4.0-2.26  OpenMPI/1.10.3 R/3.3.3
Rscript ./create_model.R $study $expr_RDS $geno $gene_annot $snp_annot \
    $n_k_folds $alpha $out_dir $chrom $snpset $window

