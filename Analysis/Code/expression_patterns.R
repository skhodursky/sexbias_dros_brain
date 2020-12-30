library(ggplot2)
library(cowplot)
library(ggsignif)
library(GenomicFeatures)
library(Biobase)
library(edgeR)
library(stringr)
library(RColorBrewer)

setwd('/Users/skhodursky/Desktop/LabWork/Sexbias_final/Analysis/')
rm(list=ls(all=TRUE)) 



####Functions for general comparison of expression levels
load_exp_data <- function(species){
   tpm <- read.csv(paste(species, '/All_Genes_TPM_', species, '_sum.csv', sep='')
                   , stringsAsFactors = FALSE)
  finite_p <- read.csv(paste(species, '/genes_with_finite_p_val.txt', sep='')
                       ,  stringsAsFactors = FALSE)
  tpm <- tpm[tpm$genes %in% finite_p$gene_id, ]
  out <- tpm
  return(out)
}

normalize <- function(tpm){
  library(edgeR)
  #normalize all using the TMM method from edgeR. 
  
  all_data <- tpm
  
  scale <- calcNormFactors(all_data[, c(2:dim(all_data)[2])])
  for(i in c(1:length(scale))){
    all_data[, (i+1)] <- all_data[, (i+1)]/scale[i]
  }
 
  all_data_norm <-cbind(data.frame(Dmel_name = all_data[, c(1)]), ((all_data[, c(2:dim(all_data)[2])])))
  return(all_data_norm)
  }

adjust_pval <- function(pval_vector){
  p <- p.adjust(c(pval_vector),method = "BH" )
  p[p <= 0.001] <- "***"
  p[0.001<p & p <= 0.01] <- "**"
  p[0.01< p & p <= 0.05] <- "*"
  p[p>0.05] <- "NS"
  return(p)
}

get_sbgenes <- function(species){
  sbg <-  read.csv(paste(species,'/sex_biased_genes_',species,'.txt', sep = '')
                   , stringsAsFactors = FALSE)
  m <- sbg[sbg$sex == 'm'  ,"gene_id"]
  f <- sbg[sbg$sex == 'f',"gene_id"]
  out <- list(m,f)
  names(out) <- c("m","f")
  return(out)
}
get_strain_names <- function(tpm){
  out <- unique(substr(colnames(tpm[,c(2:dim(tpm)[2])]),5,9))
  return(sort(out))
}
get_median_exp <- function(tpm){

 
  out <- data.frame(gene_id = tpm$Dmel_name)
 
  return(cbind(out,data.frame(exp = rowMedians(as.matrix(tpm[,c(2:dim(tpm)[2])])))))
}

comp_sbg_TPM <- function(species){

 
  all_exp <- read.csv(paste(species,'/genes_with_finite_p_val.txt',sep=''),
                      stringsAsFactors = FALSE) 
  all_exp <- all_exp$gene_id

  tpm <- load_exp_data(species)
  tpm <- normalize(tpm)
  tpm <- get_median_exp(tpm)
  ##First take median in each strain then take median of strains

  sb_genes <- get_sbgenes(species)
  m_biased <- sb_genes$m
  f_biased <- sb_genes$f
  u_biased <- all_exp[(!(all_exp %in% m_biased)) & (!(all_exp %in% f_biased)) ]
  
  file_list <- list.files("../refgtfs/")
  cond_1 <- grepl(species, file_list, ignore.case = T)
  cond_2 <- grepl("gtf", file_list, ignore.case = T)
  gtf_file <- paste("../refgtfs/", file_list[cond_1 & cond_2], sep='')
  gtf <- makeTxDbFromGFF(gtf_file)

  gen <- genes(gtf)
  gen_df <- as.data.frame(gen)
  x <- gen_df[grepl('X',gen_df$seqnames),"gene_id"]
  a <- gen_df[(grepl('2L',gen_df$seqnames)|
                 grepl('2R',gen_df$seqnames)|
                 grepl('3L',gen_df$seqnames)|
                 grepl('3R',gen_df$seqnames)|
                 (gen_df$seqnames == '4')),"gene_id"]

  m_x <- data.frame(exp = tpm[(tpm$gene_id %in% m_biased) & (tpm$gene_id %in% x), "exp"],  
                    type = "Male\nX")
  f_x <- data.frame(exp = tpm[(tpm$gene_id %in% f_biased) & (tpm$gene_id %in% x), "exp"],  
                    type = "Female\nX")
  u_x <- data.frame(exp = tpm[(tpm$gene_id %in% u_biased) & (tpm$gene_id %in% x), "exp"], 
                    type = "Unbiased\nX")

  m_a <- data.frame(exp = tpm[(tpm$gene_id %in% m_biased) & (tpm$gene_id %in% a), "exp"],  
                    type = "Male\nA")
  f_a <- data.frame(exp = tpm[(tpm$gene_id %in% f_biased) & (tpm$gene_id %in% a), "exp"],  
                    type = "Female\nA")
  u_a <- data.frame(exp = tpm[(tpm$gene_id %in% u_biased) & (tpm$gene_id %in% a), "exp"],  
                    type = "Unbiased\nA")

  p1 <- wilcox.test(m_x$exp, m_a$exp)
  p2 <- wilcox.test(f_x$exp, m_a$exp)
  p3 <- wilcox.test(u_x$exp, m_a$exp)
  p4 <- wilcox.test(f_a$exp, m_a$exp)
  p5 <- wilcox.test(u_a$exp, m_a$exp)

  p <- c(p1$p.value, p2$p.value, p3$p.value, p4$p.value, p5$p.value)
  p <- adjust_pval(p)


  plot <- ggplot(rbind(m_x, f_x, u_x, m_a, f_a, u_a), aes(x=type, y=log(exp+0.001),  fill = type) ) + theme_classic()+ 
  scale_fill_brewer(palette="Dark2") +
  geom_boxplot() + ylab('Log(TPM)')  + theme(legend.position = "none") + 
  ggtitle(species)+geom_signif(comparisons = list(c('Male\nX',  'Male\nA'), c('Female\nX',  'Male\nA'), 
                                                  c('Unbiased\nX',  'Male\nA'), c('Female\nA', 'Male\nA'), 
                                                  c('Unbiased\nA', 'Male\nA')), 
                               annotations = c(p[1], p[2], p[3], p[4], p[5]), y_position = c(10, 13, 16, 16, 13)) + xlab('') + ylim(-10, 17) 
  
return(plot)
}


##________________Functions for Comparisons of fold change and significance on different chromosomes.


####_________Basic Fold Change Analysis


get_fc_and_qval <- function(sex, region_name, species){
  sbg_df <- read.csv(paste(species,'/sex_biased_genes_',species,'.txt', sep = ''),
                     stringsAsFactors = FALSE)
  sbg_df$seqnames <- str_remove(sbg_df$seqnames,'Scf_')
  out <- data.frame(l2fc = c(),qval = c())
  
  if(region_name == 'X'){
    chroms <- 'X'}
  if(region_name == 'Autosomes'){
    chroms <- c('2L','3L','2R','3R','4')}
  
  for(chr in chroms){
    if (chr != '4'){
      #slice is the relevant rows from the sbg_df
      slice <- sbg_df[grepl(chr,sbg_df$seqnames) & (sbg_df$sex == sex),]
    
    }
    if(chr == '4'){
      slice <- sbg_df[(sbg_df$seqnames == chr) & (sbg_df$sex == sex),]
    }
    out <- rbind(out,data.frame(l2fc = abs(slice$l2fc),qval = -log10(slice$qval)))
  }
  out[,'region'] <- region_name
  return(out)
}


##Distribution log_fold change X vs A
make_plots_fcqval <- function(sex, species){
  
  data <- rbind(get_fc_and_qval(sex, 'X', species),
              get_fc_and_qval(sex, 'Autosomes', species))

  data$region <- factor(data$region)
  data$region <- factor(data$region, levels = c("X","Autosomes"))

#For the first part of the plot titles
  if(sex == 'm'){
    title_fp <- 'Male-biased'}
  if(sex == 'f'){
    title_fp <- 'Female-biased'}

#get the tops of the plots for setting ylim
  top_fc <- max(data$l2fc[is.finite(data$l2fc)])
  top_qval <- max(data$qval[is.finite(data$qval)])


  p1 <- wilcox.test(data$l2fc[data$region == "X"],data$l2fc[data$region != "X"])
  p2 <- wilcox.test(data$qval[data$region == "X"],data$qval[data$region != "X"])

  plot_fc <- ggplot(data,aes(x=region,y=l2fc, fill = (region))) +
  geom_boxplot() + ylab("Log2FC") + xlab('') + annotate("text",x=1.5,y=1.1*top_fc,label = paste("p = ",signif(p1$p.value,2))) +
  theme_classic() + scale_fill_brewer(palette="Dark2") + theme(legend.position = "none") + ggtitle(paste0(title_fp,'\nfold-change')) + ylim(0,1.2*top_fc)

  plot_qval <- ggplot(data,aes(x=region,y=qval, fill = (region))) +
  geom_boxplot() + ylab("-Log10(qval)") + xlab('') + annotate("text",x=1.5,y=1.1*top_qval,label = paste("p = ",signif(p2$p.value,2))) +
  theme_classic() + scale_fill_brewer(palette="Dark2") +theme(legend.position = "none") + ggtitle(paste0(title_fp,'\nsignificance')) + ylim(0,1.2*top_qval)

return(plot_grid(plot_fc,plot_qval,align = T))
}

