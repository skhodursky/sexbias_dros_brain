##load in data
library(car)
library(GenomicFeatures)
library(edgeR)
library(DESeq2)
library(ggplot2)
library(cowplot)
library(mgcv)


#function to load TPM data based on sex and species for oliver lab data
agg_TPM_oliver <- function(tissue, species, sex){
  #load table with sample info
  species_samples <- read.table(paste0('Oliver_lab/',species, '_samples.txt'), header = F, stringsAsFactors = F)
  colnames(species_samples) <- c('sample', 'species', 'sex','tissue','replicate')
  rel_samples <- species_samples$sample[species_samples$tissue == tissue & species_samples$sex == sex]
  n <- 0
  for(sample in rel_samples){
    n <- n + 1
    exp_data <- read.delim(paste0('Oliver_lab/', species, '_expression/',sample,'/gene_abund.tsv'), sep = '\t', stringsAsFactors = F)
    tpm_data <- exp_data[,c('Gene.ID','TPM')]
    tpm_data <- aggregate(tpm_data$TPM, by=list(Gene.ID=tpm_data$Gene.ID), FUN=sum)
    colnames(tpm_data) <- c('Gene.ID',sample)
    if(n == 1){
      total_data <- tpm_data
    }
    else{
      total_data <- merge(total_data, tpm_data, by = "Gene.ID")
    }
  }
  return(total_data)
  
  
}

agg_TPM_brain <- function(species){
  brain <- read.csv(paste0(species,'/all_genes_tpm_',species,'_sum.csv'),header = T,stringsAsFactors = F)
}

get_m_f_brain <- function(brain_data){
  
  data <- brain_data[,c(2:dim(brain_data)[2])]
  sex <- substr(colnames(data),9,9)
  male_data <- cbind(data.frame(Gene.ID = brain_data$genes),data[,sex == 'm'])
  female_data <- cbind(data.frame(Gene.ID = brain_data$genes),data[,sex == 'f'])
  out <- list(male_data,female_data)
  names(out) <- c("male","female")
  return(out)
}


#get median expression of a certain gene in a tissue in a given species
get_med_exp <- function(data){
  out <- data.frame(Gene.ID = data$Gene.ID, exp = rowMedians(as.matrix(data[, c(2:dim(data)[2])])))
  return(out)
}

#function to calculate tau
tau <- function(x){
  x <- log2(x +1) 
  x_hat <- x/max(x)
  return(sum(1-x_hat)/(length(x)-1))
}

#merge median expression data across tissues
exp_across_tis <- function(species){
  
  orthologs <- read.csv('Evolution_of_Sex_Bias/one_one_orthologs_mel_yak_sim.csv',header=T,stringsAsFactors = F)
  brain <- agg_TPM_brain(species)
  brain_m <- get_m_f_brain(brain)$male
  brain_f <- get_m_f_brain(brain)$female
 
  brain_f <- get_med_exp(brain_f[brain_f$Gene.ID %in% orthologs[,species],])
  brain_m <- get_med_exp(brain_m[brain_m$Gene.ID %in% orthologs[,species],])
  
  
  abdomen_m <- get_med_exp(agg_TPM_oliver('abdomen', species, 'male'))
  abdomen_f <- get_med_exp(agg_TPM_oliver('abdomen', species, 'female'))
  
  thorax_m <- get_med_exp(agg_TPM_oliver('thorax', species, 'male'))
  thorax_f <- get_med_exp(agg_TPM_oliver('thorax', species, 'female'))
  
  viscera_m <- get_med_exp(agg_TPM_oliver('viscera', species, 'male'))
  viscera_f <- get_med_exp(agg_TPM_oliver('viscera', species, 'female'))
  
  gonads_m <- get_med_exp(agg_TPM_oliver('gonad', species, 'male'))
  gonads_f <- get_med_exp(agg_TPM_oliver('gonad', species, 'female'))

  # merge expression across tissues together. Only keep tissues with TPM at least 1 in each tissue.
  male_t <- Reduce(function(x, y) merge(x,  y,  by = "Gene.ID"), 
                   list(brain_m,  abdomen_m,  thorax_m,  viscera_m,  gonads_m))
  male_t <- male_t[apply(male_t[, c(2:6)], 1, max)>1, ]
  female_t <- Reduce(function(x,  y) merge(x,  y,  by = "Gene.ID"),  
                     list(brain_f,  abdomen_f,  thorax_f,  viscera_f,  gonads_f))
  female_t <- female_t[apply(female_t[, c(2:6)], 1, max)>1, ]
  
  out <- list(male_t, female_t)
  names(out) <- c("male", "female")
  return(out)
}

#to convert Dyak gene IDs to Dmel gene IDs
convert_geneids <- function(yak_ids){
  orthologs <-read.csv('Evolution_of_Sex_Bias/one_one_orthologs_mel_yak_sim.csv',header=T,stringsAsFactors = F)
  mel_ids <- c()
  for(gene in yak_ids){
    mel_ids <- c(mel_ids,orthologs$Dmel[orthologs$Dyak == gene])
  }
  return(mel_ids)
}

load_tau <- function(species){
  male_tau <- read.csv(paste0('Tissue_Specificity/',species,'_male_tau.csv'),header = T, stringsAsFactors = F)
  female_tau <- read.csv(paste0('Tissue_Specificity/',species,'_female_tau.csv'),header = T, stringsAsFactors = F)
  out <- list(male_tau,female_tau)
  names(out) <- c("male_tau","female_tau")
  return(out)
}

msr_v_tau_plot <- function(data, stat, sex, species){
  if(stat == 'MSR'){
    if(sex == 'Male'){
      ylabel <- expression(atop('Log('*italic('MS'['sp']/'MS'['st']) * ')', 'male'))
  }else{
      ylabel <- expression(atop('Log('*italic('MS'['sp']/'MS'['st']) * ')', 'female'))}
  }
  if(stat == 'MS_sp'){
    if(sex == 'Male'){
      ylabel <- expression(atop('Log('*italic('MS'['sp']) * ')','male'))
    }else{
      ylabel <- expression(atop('Log('*italic('MS'['sp']) * ')','female'))
    }
  }
  if(stat == 'MS_st'){
    if(sex == 'Male'){
      ylabel <- expression(atop('Log('*italic('MS'['st']) * ')','male'))
    }else{
      ylabel <- expression(atop('Log('*italic('MS'['st']) * ')','female'))
    }
  }
  
  data[, 'Stat'] <- data[, stat]
  top <- max(log(data$Stat))
  p <- cor.test(data$tau,data[, stat], method = "spearman")
  #pval <-signif(p$p.value, 3)
  pval <- formatC(p$p.value, format = "e", digits = 3)
  #r <- signif(p$estimate, 3)
  r <- formatC(p$estimate, format = "f",digits = 3)
  if(species == 'Dmel'){
    xlabel <- expression('Tissue specificity in' ~ italic('D. mel'))
  }
  if(species == 'Dyak'){
    xlabel <- expression('Tissue specificity in' ~ italic('D. yak'))
  }
  plot <- ggplot(data,aes(x = tau, y=log(Stat))) + theme_classic() + geom_point(color = "gray", size = 0.5) + 
    annotate("text", x = 0.5, y = 1.2*top, label = paste('p = ', pval)) +
    annotate("text", x = 0.5, y = 1.5*top, label = paste('r = ', r)) + 
    ylab(ylabel) + xlab(xlabel) + 
    geom_smooth(color = "black",se=FALSE) + ylim(-10,15)
  return(plot)
  
}

###_____________________ Regress out the effects of Tau on MS ratio and variability

regress_out_tau_msr <- function(data){
 
  model<-gam(log(MSR) ~ s(tau, bs = "cs"),data=data)
  data[,'MSR'] <- residuals(model,type = "response")
  return(data)
}

regress_out_tau_msbet <- function(data){
  
  model<-gam(log(MS_sp) ~ s(tau, bs = "cs"),data=data)
  data[,'MS_sp'] <- residuals(model,type = "response")
  
  return(data)
}

regress_out_tau_mswit <- function(data){
  
  model<-gam(log(MS_st) ~ s(tau, bs = "cs"),data=data)
  data[,'MS_st'] <- residuals(model,type = "response")
  return(data)
}

######________________________________________________Differential Expression Analysis

write_list_of_samples <- function(tissue, species){
  species_samples <- read.table(paste0(species, '_samples.txt'), header = F, stringsAsFactors = F)
  colnames(species_samples) <- c('sample', 'species', 'sex','tissue','replicate')
  rel_samples <- species_samples$sample[species_samples$tissue == tissue]
  gtf_location <- paste0(species,'_expression/',rel_samples,'/',rel_samples,'.gtf')
  out <- data.frame(samples = rel_samples, gtf_loc = gtf_location)
  write.table(out,paste0(species,'_',tissue,'_samples.txt'), quote = F, row.names = F, col.names = F)
}



add_bias_fromFC <- function(sbg_df, cts, pheno_data, tissue){
  # pick a gene with highest l2fc and determine sex-bias
  # from that
  gene <- as.character(sbg_df$gene_id[order(sbg_df$l2fc, decreasing = TRUE)][1])
  l2fc_gene <- sbg_df$l2fc[sbg_df$gene_id == gene]
  # get expression of gene in males and females.
  male_exp <- rowMeans(cts[gene, pheno_data$Ids[(pheno_data$sex == 'male') & (pheno_data$tissue == tissue)]])
  female_exp <- rowMeans(cts[gene, pheno_data$Ids[(pheno_data$sex == 'female') & (pheno_data$tissue == tissue)]])
  # sex-biased data frame that contains sex-bias
  sbg_df_out <- sbg_df
  sbg_df_out$sex <- sbg_df$l2fc
  
  if(male_exp>female_exp){
    sbg_df_out[sbg_df_out$sex >0,'sex'] <- 'm'
    sbg_df_out[sbg_df_out$sex <0,'sex'] <- 'f'
  }else{
    sbg_df_out[sbg_df_out$sex >0,'sex'] <- 'f'
    sbg_df_out[sbg_df_out$sex <0,'sex'] <- 'm'
  }
  
  return(sbg_df_out)
}

run_DESeq <- function(species, tissue, signif){
  
  ###counts matrix generated by prepDE.py
  cts <- read.csv(paste0('Oliver_lab/DESeq_input/gene_count_matrix_', species,'_', tissue,'.csv'),header = TRUE)
  rownames(cts) <- cts$gene_id
  cts$gene_id <- NULL
  
  #generate phenotype data
  num_samples <- dim(cts)[2]
  pheno_data <- read.table(paste0('Oliver_lab/',species,'_samples.txt'),
                           header = FALSE, stringsAsFactors = FALSE)
  # dyak only has a single strain dmel has two strains
  if(species == "Dyak"){
    colnames(pheno_data) <- c('Ids', 'species', 'sex', 'tissue', 'rep')
    row.names(pheno_data) <- pheno_data$Ids
    pheno_data_final <- data.frame(sex = pheno_data[pheno_data$tissue == tissue, 'sex'])
    row.names(pheno_data_final) <- row.names(pheno_data[pheno_data$tissue == tissue, ])
    dds <- DESeqDataSetFromMatrix(countData = cts, colData = pheno_data_final, design = ~ sex)
  }else if(species == "Dmel"){
    colnames(pheno_data) <- c('Ids', 'species', 'sex', 'tissue', 'rep', 'strain')
    row.names(pheno_data) <- pheno_data$Ids
    pheno_data_final <- pheno_data[pheno_data$tissue == tissue, c('sex', 'strain')]
    dds <- DESeqDataSetFromMatrix(countData = cts, colData = pheno_data_final, design = ~ strain + sex)
  }
  

  
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep, ]
  dds <- DESeq(dds)
  res <- results(dds, alpha = signif)
  resOrdered <- res[order(res$pvalue), ]
  resOrdered <- as.data.frame(resOrdered)
  resOrdered$gene_id <- rownames(resOrdered)
  resOrdered <- na.omit(resOrdered)
  # rename the log2foldchange column to l2fc
  resOrdered$l2fc <- resOrdered$log2FoldChange
  resOrdered$log2FoldChange <- NULL
  
  sex_biased_genes <- resOrdered[resOrdered$padj<signif,
                                 c("gene_id", "l2fc", "padj")]
  colnames(sex_biased_genes) <- c("gene_id", "l2fc", "qval")
  # Use a sample from the species to get the unique list of genes
  sample <- colnames(cts)[1]
  sample_data <- read.delim(paste0('Oliver_lab/',species,'_expression/',sample,'/gene_abund.tsv'),sep='\t')
  sample_data <- as.data.frame(unique(sample_data[sample_data$Gene.ID ,c('Gene.ID','Reference')]))
  colnames(sample_data) <- c("gene_id","seqnames")
  
  # sbg_df will be the final dataframe that contains the sex-biased genes and relevant info
  # sbg_df <- as.data.frame(unique(sample_data[sample_data$Gene.ID ,c('Gene.ID','Reference')]))
  # colnames(sbg_df) <- c("gene_id","seqnames")
  sbg_df <- merge(sample_data, sex_biased_genes, by = "gene_id")
  sbg_df <- add_bias_fromFC(sbg_df, cts, pheno_data, tissue)
 
  write.csv(sbg_df, file = paste0('Oliver_lab/DESeq_output/sex_biased_genes_', species, '_', tissue, '.txt'),
            quote = FALSE, row.names = FALSE)
  write.csv(sample_data[sample_data$gene_id %in% rownames(resOrdered),],
            file = paste0('Oliver_lab/DESeq_output/genes_with_finite_p_val_', species, '_', tissue, '.txt'),
            quote = FALSE,row.names = FALSE)
  
}


get_obs_exp_sbgenes <- function(species, tissue){
  
  if(tissue != "brain"){
    sbg_df <- read.csv(paste0('Oliver_lab/DESeq_output/sex_biased_genes_', species, '_', tissue, '.txt'),
                       stringsAsFactors = FALSE)
    ref <- read.csv(paste0('Oliver_lab/DESeq_output/genes_with_finite_p_val_', species, '_', tissue, '.txt'),
                    stringsAsFactors = FALSE)
  }else if(tissue == "brain"){
    sbg_df <- read.csv(paste0(species, '/sex_biased_genes_', species, '.txt')
                       , stringsAsFactors = FALSE)
    ref <- read.csv(paste0(species, '/genes_with_finite_p_val.txt'),
                    stringsAsFactors = FALSE)
     }
  
    # get numbers of observed male and female biased genes
    mb_obs <- sum(sbg_df$sex == 'm')
    fb_obs <- sum(sbg_df$sex == 'f')
    # find the expressed genes (with finite p values) that are located on each chromosome
    genesX <- ref$gene_id[grepl('X', ref$seqnames)]
    genes2L <- ref$gene_id[grepl('2L', ref$seqnames)]
    genes2R <- ref$gene_id[grepl('2R', ref$seqnames)]
    genes3L <- ref$gene_id[grepl('3L', ref$seqnames)]
    genes3R <- ref$gene_id[grepl('3R', ref$seqnames)]
    genes4 <- ref$gene_id[ref$seqnames == '4']
    
    tot_A <- c(genes2L, genes2R, genes3L, genes3R, genes4)
    tot_X <- genesX
    
    frac_X <- length(tot_X)/dim(ref)[1]
    frac_A <- length(tot_A)/dim(ref)[1]
    
    # get the expected number of sex-biased genes on the X chromosome and Autosomes
    exp_mX <- mb_obs*frac_X
    exp_fX <- fb_obs*frac_X
    
    exp_mA <- mb_obs*frac_A
    exp_fA <- fb_obs*frac_A
    #find the total number of sex-biased genes that are located on each chromosome
    
    obs_mX <- sum(sbg_df$gene_id[sbg_df$sex == 'm'] %in% tot_X)
    obs_mA <- sum(sbg_df$gene_id[sbg_df$sex == 'm'] %in% tot_A)
    
    obs_fX <- sum(sbg_df$gene_id[sbg_df$sex == 'f'] %in% tot_X)
    obs_fA <- sum(sbg_df$gene_id[sbg_df$sex == 'f'] %in% tot_A)
    
    for_plot <- data.frame(Bias = c('Male-biased', 'Female-biased', 'Male-biased', 'Female-biased'),  
                           'Obs/Exp' = c(obs_mX/exp_mX, obs_fX/exp_fX, obs_mA/exp_mA, obs_fA/exp_fA), 
                      'Location' = c('X Chromosome', 'X Chromosome', 'Autosomes',  'Autosomes'), 
                      'Tissue' = tissue )
    
    for_test <- data.frame(Bias = c('Male-biased', 'Female-biased', 'Male-biased', 'Female-biased'), 
                           'Obs' = c(obs_mX, obs_fX, obs_mA, obs_fA), 
                           'Location' = c('X Chromosome', 'X Chromosome', 'Autosomes', 'Autosomes'), 
                           'Exp' = c(exp_mX, exp_fX, exp_mA, exp_fA), 
                           Tissue = tissue )
    
    return(list(for_plot, for_test))
    
  
  
}

get_sig_vs_brain <- function(brain_df, other_df){
  p_vals <- c()
  for(i in c(1:4)){
    # slice b is brain, other is other tissue
    slice_b <- brain_df[i,]
    slice_o <- other_df[i,]
    
    chisq_table <- data.frame(obs = c(slice_b$Obs,slice_o$Obs),exp = c(slice_b$Exp,slice_o$Exp))
    test <- chisq.test(chisq_table)
    p_vals <- c(p_vals,test$p.value)
  }
  # other df is the tissue being compared to brain
  return(cbind(other_df,p_vals))
}

#get ratios of observed to expected sex-biased genes on the X chromosome and autosomes
make_plot <- function(species){
  brain_ratio <-  get_obs_exp_sbgenes(species,"brain")[[1]]
  gonad_ratio <-  get_obs_exp_sbgenes(species,"gonad")[[1]]
  thorax_ratio <-  get_obs_exp_sbgenes(species,"thorax")[[1]]
  abdomen_ratio <-  get_obs_exp_sbgenes(species,"abdomen")[[1]]
  viscera_ratio <-  get_obs_exp_sbgenes(species,"viscera")[[1]]
  ratio_df <- rbind(brain_ratio, 
                    abdomen_ratio, 
                    viscera_ratio,
                    gonad_ratio, 
                    thorax_ratio)
  
  #get p values relative to brain
  
  abdomen_p <- get_sig_vs_brain(get_obs_exp_sbgenes(species,"brain")[[2]],
                                get_obs_exp_sbgenes(species,"abdomen")[[2]])
  viscera_p <- get_sig_vs_brain(get_obs_exp_sbgenes(species,"brain")[[2]],
                                  get_obs_exp_sbgenes(species,"viscera")[[2]])
  gonad_p <- get_sig_vs_brain(get_obs_exp_sbgenes(species,"brain")[[2]],
                              get_obs_exp_sbgenes(species,"gonad")[[2]])
  thorax_p <- get_sig_vs_brain(get_obs_exp_sbgenes(species,"brain")[[2]],
                               get_obs_exp_sbgenes(species,"thorax")[[2]])
  p_df <- rbind(abdomen_p,viscera_p,gonad_p,thorax_p)
  p_df[,'p_adj'] <- p.adjust(p_df$p_vals, method = 'BH')
  p_df$p_adj <- as.character(symnum(c(p_df$p_adj), corr = FALSE, 
                                    na = FALSE,cutpoints = c(0,1e-5, 1e-4, 1e-3, 1e-2, 5e-2,1),
                                    symbols = c("*****","****","***","**","*",'')))
  
  brain_p <- get_obs_exp_sbgenes(species,"brain")[[2]]
  brain_p[,'p_vals'] <- ''
  brain_p[,'p_adj'] <- ''
  p_df_final <- rbind(brain_p, p_df)
  
  plot <- ggplot(ratio_df, aes(x = Tissue,y = Obs.Exp,fill = Tissue)) + facet_grid(Location~ Bias) +
    scale_fill_brewer(palette="Dark2")+ geom_col() + theme_classic() +ylab('Observed/Expected')+ geom_text(data = p_df_final,aes( y = 1.5,label = p_adj),hjust = 0.5,vjust = 0.5,size=4)
  return(plot)
}