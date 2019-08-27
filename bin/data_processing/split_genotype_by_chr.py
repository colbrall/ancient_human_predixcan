#! /usr/bin/env python
# modified by Laura Colbran (1/3/18) to take GTEx v6p VCF input

import os
import sys
import gzip
import statistics

'''
Script to split a GTEx genotype file into multiple files by chromosome.

From commandline, first argument is the genotype file, second is the
prefix for the output files.  The suffix 'chrN.txt' will be added to the
prefix provided, where N is the chromosome number.

In splitting, script will only keep unambiguously stranded SNPs. I.e.,
no INDELs and no SNPs with polymorphisms A->T and vice-versa, or C->G
and vice-versa.

The input file is expected to be a tab-delimited text file including a
header row, where the header field is ID (for snp varID) and then a
variable number of fields with the sample ID numbers.  The first column
has the snp ID in the format of (chr_pos_refAll_effAll_build) and the
dosages are encoded on a 0-2 scale representing the number or imputed
number of the effect alleles the sample posseses.
'''

SNP_COMPLEMENT = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}

def split_genotype(geno_file, out_prefix):
    # Make output file names from prefix.
    geno_by_chr_fns = [out_prefix + '.chr' + str(i) + '.txt' for i in range(1,23)]
    # Open connection to each output file.
    geno_by_chr = [open(f, 'w') for f in geno_by_chr_fns]

    with gzip.open(geno_file, 'r') as geno:
        snps = set()
        for line in geno:
            l = line.decode("utf-8")
            if l.startswith("#"):
                if l.startswith("##"): continue #skip VCF info
                else: # Write header in each file
                    header = ["ID"]
                    for ind in range(9,len(l.split())): #pull individual IDs
                        header += ["-".join(l.split()[ind].split("-")[0:2])]
                    for f in geno_by_chr:
                        f.write('\t'.join(header) + "\n")
                    continue
            in_line = l.split()
            chr = in_line[0]
            refAllele = in_line[3]
            effectAllele = in_line[4]
            # Skip non_single letter polymorphisms
            if len(refAllele) > 1 or len(effectAllele) > 1:
                continue
            # Skip ambiguous strands
            if SNP_COMPLEMENT[refAllele] == effectAllele:
                continue
            varID = in_line[2]
            # Some snps have 2 rows for some reason. Attributes are nearly
            # identical. Only keep the first one found.
            if varID in snps:
                continue
            snps.add(varID)
            out_line = [varID]
            # pull decimal dosage
            # for ind in range(9,len(in_line)):
            #     out_line += [in_line[ind].split(':')[2]]

            # convert genotype to dosage
            for ind in range(9,len(in_line)):
               genotype = in_line[ind].split(':')[0]
               if genotype == "0/0": out_line += [0]
               elif genotype == "0/1" or genotype == "1/0": out_line += [1]
               elif genotype == "1/1": out_line += [2]
               else:
                   out_line += [-1] #if a genotype call wasn't made for that person
            # mean level imputation
            mean = statistics.mean([t for t in out_line[1:] if t != -1])
            out_line = [str(mean) if i == -1 else str(i) for i in out_line]

            # Write line to appropriate file
            index = int(chr) - 1
            geno_by_chr[index].write('\t'.join(out_line) + "\n")

    for f in geno_by_chr:
        f.close()

if __name__ == '__main__':
    genotype_file = sys.argv[1]
    out_prefix = sys.argv[2]
    split_genotype(genotype_file, out_prefix)
