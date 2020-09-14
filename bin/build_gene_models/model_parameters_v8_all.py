import os

# Define parameters----------------------------------------------------/
# Change these variables when adapting for different analyses.

# List of identifiers for each database you'll make:
# STUDY_NAMES = ['Muscle_Skeletal','Ovary']
STUDY_NAMES = ['Adipose_Subcutaneous','Adipose_Visceral_Omentum','Adrenal_Gland','Artery_Aorta','Artery_Coronary','Artery_Tibial',
        'Brain_Amygdala','Brain_Anterior_cingulate_cortex_BA24','Brain_Caudate_basal_ganglia','Brain_Cerebellar_Hemisphere','Brain_Cerebellum',
        'Brain_Cortex','Brain_Frontal_Cortex_BA9','Brain_Hippocampus','Brain_Hypothalamus','Brain_Nucleus_accumbens_basal_ganglia',
        'Brain_Putamen_basal_ganglia','Brain_Spinal_cord_cervical_c-1','Brain_Substantia_nigra','Breast_Mammary_Tissue','Cells_Cultured_fibroblasts',
        'Cells_EBV-transformed_lymphocytes','Colon_Sigmoid','Colon_Transverse','Esophagus_Gastroesophageal_Junction','Esophagus_Mucosa',
        'Esophagus_Muscularis','Heart_Atrial_Appendage','Heart_Left_Ventricle','Kidney_Cortex','Liver','Lung','Minor_Salivary_Gland',
        'Muscle_Skeletal','Nerve_Tibial','Ovary','Pancreas','Pituitary','Prostate','Skin_Not_Sun_Exposed_Suprapubic','Skin_Sun_Exposed_Lower_leg',
        'Small_Intestine_Terminal_Ileum','Spleen','Stomach','Testis','Thyroid','Uterus','Vagina','Whole_Blood']
# File names for gene and snp annotation:
GENE_ANNOTATION_FN = 'gencode.v27.annotation.gtf.gz'
SNP_ANNOTATION_FN = 'gtex_v8_vcfSNPs_hg38_refFix.txt.gz'
# List of genotype file names:
GENOTYPE_FNS = ['gtex_v8.vcf.gz']

# Model metadata/parameters. Keep all as strings:
SNPSET = 'full'
ALPHA = '0.5'
N_K_FOLDS = '10'
RSID_LABEL = 'RSID_dbSNP150'
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

# Define directories---------------------------------------------------/
INPUT_DIR = '/dors/capra_lab/data/genotype-tissue_expression_project/v8/'
INTER_DIR = '/scratch/colbrall/gtex_v8_data/'
OUTPUT_DIR = '/scratch/colbrall/v8_models_refFix/'
GENE_ANN_DIR = '/dors/capra_lab/projects/ancient_human_predixcan/data/gtex_v8_data/'
SNP_ANN_DIR = '/scratch/colbrall/'
SNP_ANN_INTER_DIR = 'snp_anno_gtexv8All_RefFix/'
GENOTYPE_INPUT_DIR = '/dors/capra_lab/data/genotype-tissue_expression_project/v8/genotypes_unpruned/'
EXPRESSION_INPUT_DIR = '/dors/capra_lab/data/genotype-tissue_expression_project/v8/expression/residual/'
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
