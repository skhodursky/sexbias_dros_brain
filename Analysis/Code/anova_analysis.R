##load in data
library(GenomicFeatures)
library(ggplot2)
library(edgeR)

setwd("~/Desktop/LabWork/Sexbias_final/Analysis/")
rm(list=ls(all=TRUE))



#load TPM of orthologs
add_dmel_ortho <- function(species, data, orthologs){
  dmel_ortho <- function(gene, species, orthologs){
    return(orthologs[orthologs[,species] == gene, 'Dmel'])
  }
  # just to be safe
  if(is.null(data$gene_id)){
    mel_gene <- as.vector(sapply(data$genes, dmel_ortho,species = species, orthologs = orthologs))
  }else{
    mel_gene <- as.vector(sapply(data$gene_id, dmel_ortho,species = species, orthologs = orthologs))
  }
  data <- cbind(data.frame(Dmel_name= mel_gene, stringsAsFactors=FALSE),data,stringsAsFactors = FALSE)
  return(data)
}

load_ortho_exp <- function(species){
  orthologs <-read.csv('Evolution_of_Sex_Bias/one_one_orthologs_mel_yak_sim.csv',header=T,stringsAsFactors = F)
  tpm <- read.csv(paste0(species, '/all_genes_tpm_', species, '_sum.csv'), stringsAsFactors = FALSE)
  tpm <- tpm[tpm$genes %in% orthologs[,species],]
  if(species == 'Dmel'){
    colnames(tpm)[colnames(tpm) == "genes"] <- "Dmel_name"
  }else{
    tpm <- add_dmel_ortho(species,tpm,orthologs)
    tpm$genes <- NULL
  }
  return(tpm)
}
#________normalize the TPM expression data.Normalizing after pulling
#out the 1:1 orthologs. Normalize the TPM using TMM
normalize <- function(dmel_tpm, dsim_tpm, dyak_tpm){
  #normalize all using the TMM method from edgeR. 
  
  all_data <- merge(dmel_tpm,dsim_tpm, by = "Dmel_name")
  all_data <- merge(all_data,dyak_tpm, by = "Dmel_name")
  scale <- calcNormFactors(all_data[,(2:dim(all_data)[2])], method = "TMM")
  all_data_norm <- all_data
  for(i in c(1:length(scale))){
    all_data_norm[,(i+1)] <- all_data[,(i+1)]/scale[i]
  }
  
  
  dmel_tpm_out <- all_data_norm[, colnames(dmel_tpm)]
  dsim_tpm_out <- all_data_norm[, colnames(dsim_tpm)]
  dyak_tpm_out <- all_data_norm[, colnames(dyak_tpm)]
  out <- list(dmel_tpm_out,dsim_tpm_out,dyak_tpm_out)
  names(out) <- c("dmel_tpm","dsim_tpm","dyak_tpm")
  return(out)
}

###process data for input into SAS (in order to run the nested random effects ANOVA analysis)
get_strains <- function(data){
  strains <- colnames(data)
  strains <- substr(strains[c(2:length(strains))],5,8)
  return(strains)
}
generate_sas_input <- function(tpm_mel, tpm_sim, tpm_yak, sex){
  genes <- c()
  for(gene in tpm_mel$Dmel_name){
    strains_mel <- get_strains(tpm_mel)
    
    strains_sim <- get_strains(tpm_sim)
   
    strains_yak <- get_strains(tpm_yak)
    
    data_mel <- data.frame(exp = as.numeric(as.vector(tpm_mel[tpm_mel$Dmel_name == gene,
                                                              c(2:dim(tpm_mel)[2])])), species = "Dmel", strains = strains_mel)
    data_sim <- data.frame(exp = as.numeric(as.vector(tpm_sim[tpm_sim$Dmel_name == gene,
                                                              c(2:dim(tpm_sim)[2])])), species = "Dsim", strains = strains_sim)
    data_yak <- data.frame(exp = as.numeric(as.vector(tpm_yak[tpm_yak$Dmel_name == gene,
                                                              c(2:dim(tpm_yak)[2])])), species = "Dyak", strains = strains_yak)

    data <- rbind(data_mel,data_yak,data_sim)
    if (max(data$exp)>1){
      # log normalize
      data$exp <- log2(data$exp + 0.001)
      # pareto normalize 
      data$exp <- (data$exp - mean(data$exp))/sqrt(sd(data$exp))
      write.table(data, paste0('SASUniversityEdition/myfolders/input_data_', sex, '/', gene,'.txt'),
                  quote = FALSE, row.names = FALSE)
      genes <- c(genes, gene)
    }
    
  }
  genes_df <- data.frame(gene = genes)
  # list of genes for reading into SAS
  write.table(genes_df, paste0('SASUniversityEdition/myfolders/input_data_',sex,'/genelist.txt'),
              quote = F, row.names = F)
}
###___________Append median TPM to ANOVA data
append_tpm <- function(anova, mel_df, sim_df, yak_df){
  tpm_mel <- c()
  tpm_sim <- c()
  tpm_yak <- c()
  for (gene in anova$gene_id){
    tpm_mel <- c(tpm_mel,median(as.numeric(as.vector(mel_df[mel_df$Dmel_name == gene,c(2:dim(mel_df)[2])]))))
    tpm_sim <- c(tpm_sim,median(as.numeric(as.vector(sim_df[sim_df$Dmel_name == gene,c(2:dim(sim_df)[2])]))))
    tpm_yak <- c(tpm_yak,median(as.numeric(as.vector(yak_df[yak_df$Dmel_name == gene,c(2:dim(yak_df)[2])]))))
  }
  anova[,'TPM_dmel'] <- tpm_mel
  anova[,'TPM_dsim'] <- tpm_sim
  anova[,'TPM_dyak'] <- tpm_yak
  return(anova)
}

##_Split male and female samples
split_malefemale <- function(input_data){
  data <- input_data[,c(2:dim(input_data)[2])]
  sex <- substr(colnames(data),9,9)
  male_data <- cbind(data.frame(Dmel_name = input_data$Dmel_name), data[, sex == 'm'])
  female_data <- cbind(data.frame(Dmel_name = input_data$Dmel_name), data[, sex == 'f'])
  out <- list(male_data,female_data)
  names(out) <- c("male","female")
  return(out)
}




####
####
####            ANALYZE and PLOT OUTPUT

library(GenomicFeatures)

load_anova_data <- function(){
  anova_mal <- read.csv('Evolution_of_Sex_Bias/Expression_Evolution/anova_male_yak_mel_sim.csv',
                        stringsAsFactors = FALSE)
  anova_fem <- read.csv('Evolution_of_Sex_Bias/Expression_Evolution/anova_female_yak_mel_sim.csv',
                        stringsAsFactors = FALSE)
  out <- list(anova_mal, anova_fem)
  names(out) <- c("anova_mal", "anova_fem")
  return(out)
}

load_sbg_data <- function(){
  orthologs <-read.csv('Evolution_of_Sex_Bias/one_one_orthologs_mel_yak_sim.csv',header=T,stringsAsFactors = F)
  
  
  sbg_mel <- read.csv('Dmel/sex_biased_genes_dmel.txt',header = T, stringsAsFactors = F)
  sbg_mel <- sbg_mel[sbg_mel$gene_id %in% orthologs$Dmel,]
  
  sbg_yak <- read.csv('Dyak/sex_biased_genes_dyak.txt', header = T, stringsAsFactors = F)
  sbg_yak <- sbg_yak[sbg_yak$gene_id %in% orthologs$Dyak,]
  sbg_yak <- add_dmel_ortho('Dyak', sbg_yak, orthologs)
  
  sbg_sim <-read.csv('Dsim/sex_biased_genes_dsim.txt',header = T, stringsAsFactors = F)
  sbg_sim <- sbg_sim[sbg_sim$gene_id %in% orthologs$Dsim,]
  sbg_sim <- add_dmel_ortho('Dsim', sbg_sim, orthologs)
  
  ref_gtf <- makeTxDbFromGFF("../refgtfs/dmel-all-r6.24.gtf")
  gen <- genes(ref_gtf)
  gen_df <- as.data.frame(gen)
  x <- gen_df[grepl('X', gen_df$seqnames), 'gene_id'] 
  a <- gen_df[grepl('2L', gen_df$seqnames) |
                grepl('2R', gen_df$seqnames) | 
                grepl('3L', gen_df$seqnames) |
                grepl('3R', gen_df$seqnames) | 
                (gen_df$seqnames == '4'), 'gene_id']
  
  m_mel <- sbg_mel[sbg_mel$sex == 'm',"gene_id"]
  f_mel <- sbg_mel[sbg_mel$sex == 'f',"gene_id"]
  
  m_sim <- sbg_sim[sbg_sim$sex == 'm',"Dmel_name"]
  f_sim <- sbg_sim[sbg_sim$sex == 'f',"Dmel_name"]
  
  m_yak <- sbg_yak[sbg_yak$sex == 'm',"Dmel_name"]
  f_yak <- sbg_yak[sbg_yak$sex == 'f',"Dmel_name"]
  
  m <- unique(c(m_mel,m_sim,m_yak))
  f <- unique(c(f_mel,f_sim,f_yak))
  switched <- m[m %in% f]
  
  # don't include switched genes in analysis
  x_m <- m[(m %in% x) & (!(m %in% switched))]
  auto_m <- m[(m %in% a) & (!(m %in% switched)) ]
  x_f <-f[(f %in% x) & (!(f %in% switched))]
  auto_f <-  f[(f %in% a) & (!(f %in% switched))]
  
  # unbiased genes aren't sex-biased in any species
  x_u <- x[(! x %in% c(m, f)) ]
  a_u <- a[(! a %in% c(m, f)) ]
  out <- list(x_m, auto_m, x_f, auto_f, x_u, a_u, switched)
  names(out) <- c("x_m", "auto_m", "x_f", "auto_f", "x_u", "auto_u", "switched")
  return(out)
}

#function to adjust p value
adjust_pval <- function(pval_vector){
  p <- p.adjust(c(pval_vector), method = "BH" )
  print(p)
  p <- as.character(symnum(p, corr = FALSE, na = FALSE,
                           cutpoints = c(0,1e-5, 1e-4, 1e-3, 1e-2, 5e-2,1),
                           symbols = c("*****","****","***","**","*",'NS')))
  return(p)
}



#function to test for the significance of MS ratio using simulation
get_pvalue_sim <- function(anova_data, geneset, numsim){
  #define sampling function
  sample_genomewide <- function(iter, anova_data, n){
    idx <- sample(c(1:dim(anova_data)[1]), n)
    return(median(anova_data$MSR[idx]))
  }
  #get number of genes sampled and the observed median ratio
  num_sampled <- sum(geneset %in% anova_data$gene_id)
  median_ratio_sampled <- median(anova_data$MSR[anova_data$gene_id %in% geneset])
  
  iteration <- c(1:numsim)
  sim_out <- sapply(iteration, sample_genomewide, anova_data = anova_data, n = num_sampled)
  
  if(median_ratio_sampled > median(anova_data$MSR)){
    p <- 2*(sum(sim_out > median_ratio_sampled)/numsim)
  }else{
    p <- 2*(sum(sim_out < median_ratio_sampled)/numsim)
  }
  return(p)
}

simsig_plots <- function(stat, plot_top, numsim = 1000000){
  library(ggplot2)
  library(ggsignif)
  library(cowplot)
  library(RColorBrewer)
  
  #load ANOVA data
  anova <- load_anova_data()
  anova_mal <- anova$anova_mal
  
  anova_fem <- anova$anova_fem
  
  #Load gene lists
  sex_biased_lists <- load_sbg_data()
  x_m <-sex_biased_lists$x_m
  x_f <-sex_biased_lists$x_f
  auto_m <-sex_biased_lists$auto_m
  auto_f <-sex_biased_lists$auto_f
  x_u <-sex_biased_lists$x_u
  auto_u <-sex_biased_lists$auto_u
  
  #get simulation p-values for each geneset
  genelists <- list(x_f, x_m, x_u, auto_f, auto_m, auto_u)
  p_mal <- c()
  p_fem <- c()
  for(geneset in genelists){
    p_mal <- c(p_mal, get_pvalue_sim(anova_mal, geneset, numsim))
    p_fem <- c(p_fem, get_pvalue_sim(anova_fem, geneset, numsim))
    
  }
  
  
  
  #multiple test adjust p-values
  pvals <- data.frame(p = adjust_pval(c(p_mal, p_fem)))
  #make data_frames for plotting
  anova_m <- rbind(data.frame(Stat = median(anova_mal[anova_mal$gene_id %in% x_f, stat]), class = 'Female X' ),
                   data.frame(Stat = median(anova_mal[anova_mal$gene_id %in% x_m, stat]), class = 'Male X' ),
                   data.frame(Stat=median(anova_mal[anova_mal$gene_id %in% x_u, stat]), class = 'Unb. X'),
                   data.frame(Stat = median(anova_mal[anova_mal$gene_id %in% auto_f, stat]), class = 'Female Auto.' ),
                   data.frame(Stat=median(anova_mal[anova_mal$gene_id %in% auto_m, stat]), class = 'Male Auto.'),
                   data.frame(Stat=median(anova_mal[anova_mal$gene_id %in% auto_u, stat]), class = 'Unb. Auto.'))
  
  anova_f <- rbind(data.frame(Stat = median(anova_fem[anova_fem$gene_id %in% x_f, stat]), class = 'Female X'),
                   data.frame(Stat = median(anova_fem[anova_fem$gene_id %in% x_m, stat]), class = 'Male X'),
                   data.frame(Stat=median(anova_fem[anova_fem$gene_id %in% x_u, stat]), class = 'Unb. X'),
                   data.frame(Stat=median(anova_fem[anova_fem$gene_id %in% auto_f, stat]), class = 'Female Auto.'),
                   data.frame(Stat=median(anova_fem[anova_fem$gene_id %in% auto_m, stat]), class = 'Male Auto.'),
                   data.frame(Stat=median(anova_fem[anova_fem$gene_id %in% auto_u, stat]), class = 'Unb. Auto.'))
  
  
  anova_m$class = factor(anova_m$class)
  anova_f$class = factor(anova_f$class)
  
  #Plot details
  ylabel <- expression('Median'~italic('MS'['sp']/'MS'['st']))
  plot_bot <- 0
  #plot_top <-  20
  
  
  y_position <- plot_top + 1
  
  ##plot for genes in males
  plot_m <- ggplot(anova_m,aes(x=class, y=(Stat), fill = class)) + geom_col()  + 
    annotate(geom="text", x=c(1:6), y=y_position, label=pvals$p[c(1:6)], size = 4) +
    geom_hline(yintercept = median(anova_mal$MSR), linetype="dashed") + 
    xlab('') + ylim(plot_bot,plot_top + 2) + ylab(ylabel)+
    ggtitle('Expression in males') +
    theme_classic()+ theme(legend.position="none")+
    scale_fill_brewer(palette="Dark2")
  
  
  
  plot_f <- ggplot(anova_f,aes(x=class, y=(Stat), fill = class)) + geom_col()  + 
    annotate(geom="text", x=c(1:6), y=y_position, label=pvals$p[c(7:12)], size = 4) +
    geom_hline(yintercept = median(anova_fem$MSR), linetype="dashed") + 
    xlab('') + ylim(plot_bot,plot_top + 2) + ylab(ylabel)+
    ggtitle('Expression in females') +
    theme_classic()+ theme(legend.position="none")+
    scale_fill_brewer(palette="Dark2")
  
    plot <- plot_grid(plot_f+theme(text = element_text(size=10.5),plot.title=element_text(size=10))+ 
                        scale_x_discrete(labels=c("Female\nX", "Male\nX", "Unbiased\nX", "Female\nA", "Male\nA", "Unbiased\nA")),
                    plot_m+theme(text = element_text(size=10.5),plot.title=element_text(size=10))+
                      scale_x_discrete(labels=c("Female\nX", "Male\nX", "Unbiased\nX", "Female\nA", "Male\nA", "Unbiased\nA")),
                    nrow = 1,align = T)
  return(plot)
}

#Correlate MS ratio in males and females
correlation_plots <- function(){
  library(ggplot2)
  library(ggsignif)
  library(cowplot)
  anova <- load_anova_data()
  anova_mal <- anova$anova_mal
  anova_fem <- anova$anova_fem
  sex_biased_lists <- load_sbg_data()
  x_m <-sex_biased_lists$x_m
  x_f <-sex_biased_lists$x_f
  auto_m <-sex_biased_lists$auto_m
  auto_f <-sex_biased_lists$auto_f
  x_u <-sex_biased_lists$x_u
  auto_u <-sex_biased_lists$auto_u
  
  tot_anova <- merge(anova_mal[,c('gene_id','MSR')],anova_fem[,c('gene_id','MSR')],by='gene_id')
  colnames(tot_anova) <- c("gene_id",'MSR_m','MSR_f')
  
  x_m_df <- cbind(tot_anova[tot_anova$gene_id %in% x_m,], data.frame(type = 'Male X'))
  x_f_df <- cbind(tot_anova[tot_anova$gene_id %in% x_f,], data.frame(type = 'Female X'))
  x_u_df <- cbind(tot_anova[tot_anova$gene_id %in% x_u,], data.frame(type = 'Unbiased X'))
  
  auto_m_df <- cbind(tot_anova[tot_anova$gene_id %in% auto_m,], data.frame(type = 'Male Autosomes'))
  auto_f_df <- cbind(tot_anova[tot_anova$gene_id %in% auto_f,], data.frame(type = 'Female Autosomes'))
  auto_u_df <- cbind(tot_anova[tot_anova$gene_id %in% auto_u,], data.frame(type = 'Unbiased Autosomes'))
  
  main_df <- rbind(x_m_df, x_f_df, x_u_df, auto_m_df, auto_f_df, auto_u_df)
  
  x_m_cor <- data.frame(cor = cor(log(x_m_df$MSR_m),log(x_m_df$MSR_f), method='spearman'), type ='Male X' )
  x_f_cor <- data.frame(cor = cor(log(x_f_df$MSR_m),log(x_f_df$MSR_f), method='spearman'), type ='Female X' )
  x_u_cor <- data.frame(cor = cor(log(x_u_df$MSR_m),log(x_u_df$MSR_f), method='spearman'), type ='Unbiased X' )
  
  a_m_cor <- data.frame(cor = cor(log(auto_m_df$MSR_m),log(auto_m_df$MSR_f), method='spearman'), type ='Male Autosomes' )
  a_f_cor <- data.frame(cor = cor(log(auto_f_df$MSR_m),log(auto_f_df$MSR_f), method='spearman'), type ='Female Autosomes' )
  a_u_cor <- data.frame(cor = cor(log(auto_u_df$MSR_m),log(auto_u_df$MSR_f), method='spearman'), type ='Unbiased Autosomes' )
  
  cor_df <- rbind(x_f_cor,x_m_cor,x_u_cor,a_f_cor,a_m_cor,a_u_cor)
  ylabel <- expression('Log('*italic('MS'['sp']/'MS'['st']) * '), female')
  xlabel <- expression('Log('*italic('MS'['sp']/'MS'['st']) * '), male')
  
  #show three digits behind decimal even if las digit is 0
  cor_df$cor <- paste0('r = ', format(round(cor_df$cor,3), nsmall = 3))
  main_df$type <- factor(main_df$type, levels = c("Female X","Male X","Unbiased X","Female Autosomes","Male Autosomes","Unbiased Autosomes"))
  plot <- ggplot(main_df,aes(x=log(MSR_m),y=log(MSR_f))) + geom_point() + geom_text(data =cor_df,aes(x = -3, y = 8,label = cor),hjust = 0.5,vjust = 0.5,size=3)+ 
    facet_wrap(~type) +theme_classic() + xlab(xlabel) + ylab(ylabel)+ coord_fixed()
  return(plot)
}


