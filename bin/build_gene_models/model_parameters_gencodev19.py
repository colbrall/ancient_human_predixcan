import os

# Define parameters----------------------------------------------------/
# Change these variables when adapting for different analyses.

# List of identifiers for each database you'll make:
STUDY_NAMES = ['Muscle_Skeletal','Ovary','Whole_Blood']
# File names for gene and snp annotation:
GENE_ANNOTATION_FN = 'gencode.v19.genes.patched_contigs.gtf.gz'
SNP_ANNOTATION_FN = 'GTEx_OMNI_genot_1KG_imputed_var_info4_maf01_CR95_CHR_POSb37_ID_REF_ALT_release_v6.txt.gz'
# List of genotype/expression file names:
GENOTYPE_FNS = ['GTEx_Analysis_20150112_OMNI_2.5M_5M_450Indiv_chr1to22_genot_imput_info04_maf01_HWEp1E6_ConstrVarIDs.vcf.gz']
EXPRESSION_RPKM = ['GTEx_Data_20150112_RNAseq_RNASeQCv1.1.8_gene_rpkm.gct.gz']
EXPRESSION_READS = ['GTEx_Data_20150112_RNAseq_RNASeQCv1.1.8_gene_reads.gct.gz']
SAMPLE_ATTRIBUTES = ['/dors/capra_lab/data/genotype-tissue_expression_project/v6p/phs000424.v6.pht002743.v6.p1.c1.GTEx_Sample_Attributes.GRU.txt.gz']
COVARIATES = "/dors/capra_lab/data/genotype-tissue_expression_project/v6p/GTEx_Analysis_v6p_eQTL_covariates/"

# Model metadata/parameters. Keep all as strings:
SNPSET = '1kG'
ALPHA = '0.5'
N_K_FOLDS = '10'
RSID_LABEL = 'RSID_dbSNP137'
WINDOW = '1e6'

# Names for intermediate files-----------------------------------------/
# File name of output for parse_gtf.py:
GENE_ANNOT_INTER1 = GENE_ANNOTATION_FN[:-6] + 'parsed.txt'
# File name of output for geno_annot_to_RDS.R:
GENE_ANNOT_INTER2 = GENE_ANNOT_INTER1[:-3] + 'RDS'
# File name prefix of outputs from split_snp_annot_by_chr.py:
SNP_ANN_INTER_PREFIX1 = SNP_ANNOTATION_FN[:-7]
# File name prefix for input files to snp_annot_to_RDS.R:
SNP_ANN_INTER_PREFIX2 = SNP_ANN_INTER_PREFIX1 + '.chr'
# File name prefixes for output files from split_genotype_by_chr.py:
GENOTYPE_INTER_PREFIX = map(lambda x: x[:-7], GENOTYPE_FNS)
# File names for output files from expr_to_transposed_RDS.R:
#EXPR_INTER = map(lambda x: x[:-3] + "RDS", EXPRESSION_RPKM)

# Define directories---------------------------------------------------/
INPUT_DIR = '../../'
INTER_DIR = '../../data/gtex_data/'
OUTPUT_DIR = '../../data/models/'
GENE_ANN_DIR = 'phe000006.v2.GTEx_RNAseq.marker-info.MULTI/'
SNP_ANN_DIR = '/dors/capra_lab/data/predixcan_models/gtex_v6_1kg/'
SNP_ANN_INTER_DIR = 'snp_ann/'
GENOTYPE_INPUT_DIR = 'phg000520.v2.GTEx_MidPoint_Imputation.genotype-calls-vcf.c1/'
EXPRESSION_INPUT_DIR = 'phe000006.v2.GTEx_RNAseq.expression-data-matrixfmt.c1/'
GENOTYPE_INTER_DIR = 'genotypes/'
EXPRESSION_INTER_DIR = 'expression/'
MODEL_BY_CHR_DIR = INTER_DIR + 'model_by_chr/'
HOME_DIR = os.path.dirname(os.path.realpath(__file__))
ALL_BETAS_FILES = list(map(lambda x: OUTPUT_DIR + 'allBetas/' + x + '.allBetas.txt', STUDY_NAMES))
ALL_COVARIANCES_FILES = list(map(lambda x: OUTPUT_DIR + 'allCovariances/' + x + '_' + SNPSET + '_alpha' + ALPHA + '_window' + WINDOW + '.txt', STUDY_NAMES))
ALL_LOGS_FILES = list(map(lambda x: OUTPUT_DIR + 'allLogs/' + x + '.allLogs.txt', STUDY_NAMES))
ALL_META_DATA_FILES = list(map(lambda x: OUTPUT_DIR + 'allMetaData/' + x + '.allMetaData.txt', STUDY_NAMES))
ALL_RESULTS_FILES = list(map(lambda x: OUTPUT_DIR + 'allResults/' + x + '.allResults.txt', STUDY_NAMES))
DB_FILES = list(map(lambda x: OUTPUT_DIR + 'dbs/' + x + '_' + SNPSET + '_alpha' + ALPHA + '_window' + WINDOW + '.db', STUDY_NAMES))
FILTERED_DB_FILES = list(map(lambda x: x[:-3] + '_filtered.db', DB_FILES))
