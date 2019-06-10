# reads normalized expression bed files downloaded from GTEx portal, reformats,
# then saves as RDS

library('readr')

makeOutput <- function(tiss_name, norm_path, cov_path, out_dir) {
    #print("writing output...")
    norm <- read.table(norm_path, stringsAsFactors = FALSE,
        header = TRUE)

    # set rownames to Gene ids, remove non-expression data columns
    rownames(norm) <- norm[,4]
    norm[,-c(1,2,3,4)]

     # truncate ind IDs for matching to genotype
    ids <- row.names(norm)
    for (i in 1:length(ids)) {
        ids[i] <- paste(strsplit(ids[i],'[.]')[[1]][1:2],collapse = "-")
    }
    row.names(norm) <- ids

    # Correct for covariates if they were provided
    if (!is.na(cov_path)) {
        covariate <- read.table(paste(cov_path,cov_files[i],sep=""), stringsAsFactors = FALSE,
                                  header = TRUE, row.names = 1)
        covariate <- t(covariate)

        # filter to set of people in Eric's covariate file (accounts for some filtering I couldn't replicate independently)
        ids <- row.names(covariate)
        for (i in 1:length(ids)) {
            ids[i] <- paste(strsplit(ids[i],'[.]')[[1]][1:2],collapse = "-")
        }
        row.names(covariate) <- ids
        norm <- norm[rownames(norm) %in% rownames(covariate),]

        # confirm same order
        norm[ order(row.names(norm)),]
        covariate[ order(row.names(covariate)),]

        for (i in 1:length(colnames(norm))) {
            fit <- lm(norm[,i] ~ covariate) #gene ~ covariates
            norm[,i] <- fit$residuals
        }
    }
    saveRDS(norm, paste0(paste0(out_dir,tiss_name,""),".RDS",""))
    write.table(t(norm), file=paste0(paste0(out_dir,tiss_name,""),".txt",""),quote=F,col.names=T,row.names=T,sep='\t') #for the metadata script
}

argv <- commandArgs(trailingOnly = TRUE)
tiss<- argv[1]
norm_file <- paste(argv[2],"Analysis.v6p.normalized.expression.bed.gz",sep="_")
out_dir <- argv[3]

# handle possible presence of covariate file
covariate_file <- ifelse(length(argv) == 4, paste(argv[4],"Analysis.v6p.covariates.txt",sep="_"), NA)

print(tiss)
makeOutput(tiss,norm_file,covariate_file,out_dir)
