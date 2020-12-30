#####Load in TPM or FPKM data into DF.
library(stringr)

rm(list=ls(all=TRUE))

# function to load expression data
load_exp_data <- function(directory){
  list <- list.files(directory)
  n <- 0
  for (s in list){
    n <- n + 1
    print(n)
    data <- read.delim(paste0(directory, '/', s, '/gene_abund.txt'),
                       header = TRUE, stringsAsFactors = F)
    # if there are multiple instances of a given gene sum TPM across
    tpm_aggregated <- aggregate(data$TPM, by=list(Gene.ID=data$Gene.ID), FUN=sum)
    colnames(tpm_aggregated) <- c("genes", s)
    if (n == 1){
      tpm_all <- tpm_aggregated
    }
    else{
      tpm_all <- merge(tpm_all, tpm_aggregated, by = "genes")
    }
   }
  return(tpm_all)
}
