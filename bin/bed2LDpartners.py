#!/usr/bin/env python

"""
bed2LDpartners.py - David Rinker 2019

 Change log:
5/31/19 Laura Colbran-- updated fro python3
 -----------------------------------------------------------------------------

"""

usage = """
USAGE:
bed2LDpartners.py POP RSQUARED SNP.bed

POP: 1KG POPULATION ("EAS", "EUR", "SAS")
RSQUARED: minimum LD (0.5 - 1.0)
SNP.bed: path/to/querySNP.bed

"""

import os
import sys
import pymysql

POP = sys.argv[1]
# POP = "EUR"

minLDrequirment = sys.argv[2]
# minLDrequirment= 0.80

###_Input SNP list as a bed file

SNP_FILE = sys.argv[3]
# SNP_FILE = "test.snp.bed"

###FUNCTIONS####

def bed_2_ChrPos (in_file):
	out_list=[]
	for line in open(in_file):
		if not line.startswith("#"):
			seg = line.rstrip().split('\t')
			CHR = seg[0].split('r')[1]
			POS = seg[2]
			out_list += [(CHR,POS)]

	return out_list

def msqlquery (MYSQLTABLE, RSQ, POS):
	ldSNPs = []
	cursor = cnx.cursor()

	query = ("SELECT POS1 FROM "+ MYSQLTABLE + " WHERE (CHROM = %s) AND (Rsquared >= %s) AND (POS2 = %s)")
	cursor.execute(query, (CHR,RSQ,POS))
	for row in cursor:
		ldSNPs.append(str(row['POS1']))

	query = ("SELECT POS2 FROM "+ MYSQLTABLE + " WHERE (CHROM = %s) AND (Rsquared >= %s) AND (POS1 = %s)")
	cursor.execute(query, (CHR,RSQ,POS))
	for row in cursor:
		ldSNPs.append(str(row['POS2']))

	return ldSNPs

#################################
##_Get LD partners for SNPs in Input SNP from pre-computed 1000 Genomes LD database
##_IMPORTANT NOTE: The LD database was constructed using +/- 50kb windows around queried SNP and only SNPs with Rsquared > 0.5 were retained
##for reference, vcftools command line was: vcftools --gzvcf AFR.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --hap-r2 --min-r2 0.5 --ld-window-bp 500000 --out AFR.chr1.phase3.ld0.5.collection

if POP=='EUR':
	DATABASETABLE='SNP_LD_1KG_PHASE3v5'
else:
	DATABASETABLE='SNP_LD_1KG_PHASE3v5_' + POP

snp_input_list = bed_2_ChrPos(SNP_FILE)

SNPout = []

cnx = pymysql.connect(host='vgi01.accre.vanderbilt.edu',
							 user='rinkerd',
							 password='wretched-calculator',
							 db='1kg_snpld_phase3',
							 charset='utf8mb4',
							 cursorclass=pymysql.cursors.DictCursor)
for snp in snp_input_list:
	CHR=snp[0]
	POS=int(snp[1])
	ldsnps = []

	ldsnps = msqlquery(DATABASETABLE, minLDrequirment, POS)

	# print '[%s]' % ', '.join(map(str, ldsnps))

	if not ldsnps:
		continue
	else:
		for ldsnp in ldsnps:
			SNPout += [(CHR,str(ldsnp))]

cnx.close()

##CONSOLIDATE OUTPUT AND PRINT IN BED FORMAT
# print(set(SNPout + snp_input_list))
uniqueSNPlist = sorted(set(SNPout + snp_input_list))

for t in uniqueSNPlist:
	print('chr%s\t%s\t%s' % (t[0], int(t[1])-1, t[1]))
