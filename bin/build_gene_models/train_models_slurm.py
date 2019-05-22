#!/usr/bin/env python
# @author Laura Colbran, 12/2018
# adapted from train_models.py to run on Slurm with Python 3

import subprocess
import time

from model_parameters import *

CMD = 'sbatch --export=study={0},expr_RDS={1},geno={2},gene_annot={3},snp_annot={4},' + \
    'n_k_folds={5},alpha={6},out_dir={7},chrom={8},snpset={9},window={10} ' + \
    '-J {0}_model_chr{8} -D {11} --output=%x.%j.out train_model_by_chr.slurm'

geno_prefixes = list(GENOTYPE_INTER_PREFIX)
gene_annot = INTER_DIR + GENE_ANNOT_INTER2
for i, study in enumerate(STUDY_NAMES):
    expression_RDS = INTER_DIR + EXPRESSION_INTER_DIR + study + ".RDS"
    #print(i,study)
    for chr in range(1,23):
        #print(str(chr))            
        geno = INTER_DIR + GENOTYPE_INTER_DIR + geno_prefixes[0] + '.chr' + str(chr) + '.txt'
        snp_annot = INTER_DIR + SNP_ANN_INTER_DIR + SNP_ANN_INTER_PREFIX2 + str(chr) + '.RDS'
        cmd = CMD.format(study,expression_RDS,geno,gene_annot,snp_annot,
            N_K_FOLDS,ALPHA,MODEL_BY_CHR_DIR,str(chr),SNPSET,WINDOW,HOME_DIR)
        print(cmd)
        subprocess.call(cmd, shell=True)
        time.sleep(2)

