##load in data
library(car)
library(GenomicFeatures)
library(edgeR)
rm(list=ls(all=TRUE))
setwd('/Users/skhodursky/Desktop/LabWork/Sexbias_final/Sex_Biased_Genes_DESeq/')

###Load Expression Data

load_exp_data <- function(species){
brain <- read.csv(paste0(species,'/all_genes_tpm_',species,'_sum.csv'),header = T,stringsAsFactors = F)

abdomen <- read.csv(paste0('~/Desktop/LabWork/OliverLab/',species,'/all_genes_tpm_',species,'_Abdomen.csv'),header = T,stringsAsFactors = F)
abdomen <- subset(abdomen,select = -c(X))
thorax <- read.csv(paste0('~/Desktop/LabWork/OliverLab/',species,'/all_genes_tpm_',species,'_Thorax.csv'),header = T,stringsAsFactors = F)
thorax <- subset(thorax,select = -c(X))
digestive<- read.csv(paste0('~/Desktop/LabWork/OliverLab/',species,'/all_genes_tpm_',species,'_Digestive.csv'),header=T,stringsAsFactors = F)
digestive <- subset(digestive,select = -c(X))
gonads <- read.csv(paste0('~/Desktop/LabWork/OliverLab/',species,'/all_genes_tpm_',species,'_Gonads.csv'),header=T,stringsAsFactors = F)
gonads <- subset(gonads,select = -c(X))

out <- list(brain,abdomen,thorax,digestive,gonads)
names(out) <- c("brain","abdomen","thorax","digestive","gonads")
return(out)
}

get_malefemale_brain <- function(brain_data){
  data <- brain_data[,c(2:dim(brain_data)[2])]
  sex <- substr(colnames(data),9,9)
  male_data <- cbind(data.frame(genes = brain_data$genes),data[,sex == 'm'])
  female_data <- cbind(data.frame(genes = brain_data$genes),data[,sex == 'f'])
  out <- list(male_data,female_data)
  names(out) <- c("male","female")
  return(out)
}
get_malefemale_oliver <- function(oliver_data){
  data <- oliver_data[,c(2:dim(oliver_data)[2])]
  sex <- substr(colnames(data),12,12)
  male_data <- cbind(data.frame(genes = oliver_data$genes),data[,sex == 'm'])
  female_data <- cbind(data.frame(genes = oliver_data$genes),data[,sex == 'f'])
  out <- list(male_data,female_data)
  names(out) <- c("male","female")
  return(out)
  
}
get_median_exp <- function(data){
  out <- data.frame(gene_id = data$genes,exp = rowMedians(as.matrix(data[,c(2:dim(data)[2])])))
  return(out)
}

compile_expression_data <- function(species){
exp <- load_exp_data(species)
brain <- exp$brain
abdomen <- exp$abdomen
thorax <- exp$thorax
gonads <- exp$gonads
digestive <- exp$digestive
orthologs <- read.csv('Evolution_of_Sex_Bias/one_one_orthologs_mel_yak_sim.csv',header=T,stringsAsFactors = F)
brain <- get_malefemale_brain(brain)
brain_f <- get_median_exp(brain$female)
brain_m <- get_median_exp(brain$male)
brain_f <- brain_f[brain_f$gene_id %in% orthologs[,species],]
brain_m <- brain_m[brain_m$gene_id %in% orthologs[,species],]


abdomen <- get_malefemale_oliver(abdomen)
abdomen_f <- get_median_exp(abdomen$female)
abdomen_m <- get_median_exp(abdomen$male)

thorax <- get_malefemale_oliver(thorax)
thorax_f <- get_median_exp(thorax$female)
thorax_m <- get_median_exp(thorax$male)

digestive <- get_malefemale_oliver(digestive)
digestive_f <- get_median_exp(digestive$female)
digestive_m <- get_median_exp(digestive$male)

gonads <- get_malefemale_oliver(gonads)
gonads_f <- get_median_exp(gonads$female)
gonads_m <- get_median_exp(gonads$male)

male_t <- Reduce(function(x, y) merge(x, y, by = "gene_id"), list(brain_m, abdomen_m, thorax_m,digestive_m,gonads_m))
male_t <- male_t[apply(male_t[,c(2:6)],1,max)>1,]
female_t <- Reduce(function(x, y) merge(x, y, by = "gene_id"), list(brain_f, abdomen_f, thorax_f,digestive_f,gonads_f))
female_t <- female_t[apply(female_t[,c(2:6)],1,max)>1,]

out <- list(male_t,female_t)
names(out) <- c("male","female")
return(out)
}

tau <- function(x){
  x <- log2(x +1) 
  x_hat <- x/max(x)
  sum(1-x_hat)/(length(x)-1)
}

brain_highest <- function(x){
  which.max(x) == 1
}
convert_geneids <- function(yak_ids){
  orthologs <-read.csv('Evolution_of_Sex_Bias/one_one_orthologs_mel_yak_sim.csv',header=T,stringsAsFactors = F)
  mel_ids <- c()
  for(gene in yak_ids){
    mel_ids <- c(mel_ids,orthologs$Dmel[orthologs$Dyak == gene])
  }
 return(mel_ids)
}

load_tau_by_species <- function(species){
  male_tau <- read.csv(paste0('Tissue_Specificity/',species,'_male_tau.csv'),header = T, stringsAsFactors = F)
  female_tau <- read.csv(paste0('Tissue_Specificity/',species,'_female_tau.csv'),header = T, stringsAsFactors = F)
  out <- list(male_tau,female_tau)
  names(out) <- c("male_tau","female_tau")
  return(out)
}
####Functions for plotting F vs tau, or ms_bet vs tau, or ms_wit vs tau
f_v_tau_subplot <- function(data,stat,sex, species){
  if(stat == 'MS_sp'){type <- 'Log(Inter. var.)'}
  if(stat == 'MS_st'){type <- 'Log(Intra. var.)'}
  if(stat == 'f_val'){type <- 'Log(F-statistic)'}
  
  data[,'Stat'] <- data[,stat]
  top <- max(log(data$Stat))
  p <- cor.test(data$tau,data[,stat],method = "spearman")
  pval <-signif(p$p.value,3)
  r <- signif(p$estimate,3)
  if(sex == 'Male'){
  ylabel <- expression('Log('~italic('MS'['sp']/'MS'['st']) ~ '), male')
  }
  if(sex == 'Female'){
    ylabel <- expression('Log('~italic('MS'['sp']/'MS'['st']) ~ '), female')
  }
  if(species == 'Dmel'){
    xlabel <- expression('Tissue specificity in' ~ italic('D. mel'))
  }
  if(species == 'Dyak'){
    xlabel <- expression('Tissue specificity in' ~ italic('D. yak'))
  }
  plot <- ggplot(data,aes(x = tau,y=log(Stat))) + theme_classic() + geom_point(color = "gray",size = 0.5) + 
    annotate("text", x = 0.5, y = 1.2*top, label = paste('p = ',pval)) +
    annotate("text", x = 0.5, y = 1.5*top, label = paste('r = ',r)) + 
    ylab(ylabel) + xlab(xlabel) + 
  geom_smooth(color = "black",se=FALSE)
  return(plot)
  
}

###_____________________ Regress out the effects of Tau on F and variability
load_sbg_data <- function(){
  orthologs <-read.csv('Evolution_of_Sex_Bias/one_one_orthologs_mel_yak_sim.csv',header=T,stringsAsFactors = F)
  
  
  sbg_mel <- read.csv('Dmel/sex_biased_genes_dmel_includingStrain.txt',header = T, stringsAsFactors = F)
  sbg_mel <- sbg_mel[sbg_mel$gene_id %in% orthologs$Dmel,]
  sbg_yak <- read.csv('Dyak/sex_biased_genes_dyak_includingStrain.txt', header = T, stringsAsFactors = F)
  sbg_yak <- sbg_yak[sbg_yak$gene_id %in% orthologs$Dyak,]
  sbg_sim <-read.csv('Dsim/sex_biased_genes_dsim_includingStrain.txt',header = T, stringsAsFactors = F)
  sbg_sim <- sbg_sim[sbg_sim$gene_id %in% orthologs$Dsim,]
  dmel_g <- c()
  for (gene in sbg_yak$gene_id){
    
    dmel_g <- c(dmel_g,orthologs[orthologs$Dyak == gene,'Dmel'])
  }
  sbg_yak[,"Gene_dmel"] <- dmel_g
  dmel_g <- c()
  for (gene in sbg_sim$gene_id){
    
    dmel_g <- c(dmel_g,orthologs[orthologs$Dsim == gene,'Dmel'])
  }
  sbg_sim[,"Gene_dmel"] <- dmel_g
  
  ref_gtf <- makeTxDbFromGFF("../refgtfs/dmel-all-r6.24.gtf")
  gen <- genes(ref_gtf)
  gen_df <- as.data.frame(gen)
  x <- gen_df[grepl('X',gen_df$seqnames),'gene_id'] 
  a <- gen_df[grepl('2L',gen_df$seqnames) | grepl('2R',gen_df$seqnames) | grepl('3L',gen_df$seqnames) |grepl('3R',gen_df$seqnames) | gen_df$seqnames == '4','gene_id']
  
  m_mel <- sbg_mel[sbg_mel$sex == 'm',"gene_id"]
  f_mel <- sbg_mel[sbg_mel$sex == 'f',"gene_id"]
  
  m_sim <- sbg_sim[sbg_sim$sex == 'm',"Gene_dmel"]
  f_sim <- sbg_sim[sbg_sim$sex == 'f',"Gene_dmel"]
  
  m_yak <- sbg_yak[sbg_yak$sex == 'm',"Gene_dmel"]
  f_yak <- sbg_yak[sbg_yak$sex == 'f',"Gene_dmel"]
  
  m <- unique(c(m_mel,m_sim,m_yak))
  f <- unique(c(f_mel,f_sim,f_yak))
  switched <- m[m %in% f]
  
  x_m <- m[(m %in% x) & (!(m %in% switched))]
  auto_m <- m[(m %in% a) & (!(m %in% switched)) ]
  x_f <-f[(f %in% x) & (!(f %in% switched))]
  auto_f <-  f[(f %in% a) & (!(f %in% switched))]
  
  
  x_u <- x[(! x %in% c(m,f)) ]
  a_u <- a[(! a %in% c(m,f)) ]
  out <- list(x_m,auto_m,x_f,auto_f,x_u,a_u,switched)
  names(out) <- c("x_m","auto_m","x_f","auto_f","x_u","auto_u","switched")
  return(out)
}



compare_gene_classes <- function(anova,stat,class1,class2){
  p <- wilcox.test(anova[anova$gene_id %in% class1,stat],anova[anova$gene_id %in% class2,stat])
  return(p$p.value)
}
adjust_pval <- function(pval_vector){
  p <- p.adjust(c(pval_vector),method = "BH" )
  #p <- pval_vector
  p <- as.character(symnum(p, corr = FALSE, na = FALSE,cutpoints = c(0,1e-5, 1e-4, 1e-3, 1e-2, 5e-2,1),symbols = c("*****","****","***","**","*",'NS')))
  return(p)
}
se <- function(data){
  return(sd(data)/sqrt(length(data)))
}





regress_out_tau_fval <- function(data){
library(mgcv)
model<-gam(log(f_val) ~ s(tau, bs = "cs"),data=data)
data[,'f_val'] <- residuals(model,type = "response")
return(data)
}

regress_out_tau_msbet <- function(data){
  library(mgcv)
  
  model<-gam(log(MS_sp) ~ s(tau, bs = "cs"),data=data)
  data[,'MS_sp'] <- residuals(model,type = "response")
  
  return(data)
}

regress_out_tau_mswit <- function(data){
  library(mgcv)
  model<-gam(log(MS_st) ~ s(tau, bs = "cs"),data=data)
  data[,'MS_st'] <- residuals(model,type = "response")
  return(data)
}


