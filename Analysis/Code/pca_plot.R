library(ggfortify)
rm(list=ls(all=TRUE))
setwd('/Users/skhodursky/Desktop/LabWork/Sexbias_final/Sex_Biased_Genes_DESeq/')

data <- read.csv('Dsim/all_genes_tpm_dsim_sum.csv', header = T, stringsAsFactors = F)
autoplot(prcomp(t(data[,c(2:11)])), label = T, label.size = 2) + theme_classic()