##Look at the direction of selection
library(phytools)
library(ape)
library(edgeR)
library(GenomicFeatures)

rm(list=ls(all=TRUE)) 
setwd('/Users/skhodursky/Desktop/LabWork/Sexbias_final/Analysis/')

# load anova data
load_anova_data <- function(){
  anova_mal <- read.csv('Evolution_of_Sex_Bias/Expression_Evolution/anova_male_yak_mel_sim.csv', stringsAsFactors = FALSE)
  anova_fem <- read.csv('Evolution_of_Sex_Bias/Expression_Evolution/anova_female_yak_mel_sim.csv', stringsAsFactors = FALSE)
  out <- merge(anova_mal, anova_fem, by = "gene_id")
  
  return(out)
}

get_sbgenes <- function(species, x_or_a){
  sbg <-  read.csv(paste(species, '/sex_biased_genes_', species, '.txt', sep = ''), stringsAsFactors = FALSE)
  if(x_or_a == 'X'){
    m <- sbg[sbg$sex == 'm'  & grepl('X',sbg$seqnames), "gene_id"]
    f <- sbg[sbg$sex == 'f'  & grepl('X',sbg$seqnames), "gene_id"]
  }else{
    m <- sbg[sbg$sex == 'm'  & (grepl('2L',sbg$seqnames)|
                                  grepl('2R',sbg$seqnames)|
                                  grepl('3L',sbg$seqnames)|
                                  grepl('3R',sbg$seqnames)|
                                  (sbg$seqnames == '4')), "gene_id"]
    f <- sbg[sbg$sex == 'f'  & (grepl('2L',sbg$seqnames)|
                                  grepl('2R',sbg$seqnames)|
                                  grepl('3L',sbg$seqnames)|
                                  grepl('3R',sbg$seqnames)|
                                  sbg$seqnames == '4'), "gene_id"]
    
  }
  
  
  ortho <- read.csv("Evolution_of_Sex_Bias/one_one_orthologs_mel_yak_sim.csv", stringsAsFactors = FALSE)
  colnames(ortho) <-c("Dmel", "Dsim", "Dyak")
  #convert to Dmel names
  m <- ortho[ortho[,species] %in% m, "Dmel"]
  f <- ortho[ortho[,species] %in% f, "Dmel"]
  #Make sure that the gene is in the anova matrix
  anova = load_anova_data()
  m <- m[m %in% anova$gene_id]
  f <- f[f %in% anova$gene_id]
  out <- list(m, f)
  names(out) <- c("m", "f")
  return(out)
}

get_direction_of_gexp_change <- function(sex, species, anova , x_or_a){
  fc_sb_m <- c()
  fc_sb_f <- c()
  lists <-get_sbgenes(species, x_or_a)
  sex_biased <- unlist(lists[sex])
  
  for(gene in sex_biased){
    
    if(gene %in% anova$gene_id){
      
      m <- as.numeric(anova[anova$gene_id == gene,c("TPM_dmel.x", "TPM_dsim.x", "TPM_dyak.x")])
      m <- log2(m + 1)
      names(m) <- c("Dmel", "Dsim", "Dyak")
      
      f <- as.numeric(anova[anova$gene_id == gene,c("TPM_dmel.y", "TPM_dsim.y", "TPM_dyak.y")])
      f <- log2(f + 1)
      names(f) <- c("Dmel", "Dsim", "Dyak")
      anc_m <- get_ancState(m)
      
      anc_f <- get_ancState(f)
      m1 <- (m[species]-abs(anc_m))/abs(anc_m)
      f1 <- (f[species]-abs(anc_f))/abs(anc_f)
      fc_sb_m<-c(fc_sb_m, m1)
      fc_sb_f<-c(fc_sb_f, f1)
      
    }
  }
  return(list(fc_sb_m, fc_sb_f))
}
get_ancState <- function(vec){
  tree<-read.tree(file='Evolution_of_Sex_Bias/Drosophila_Tree/drosophila-25spec-time-tree.tre')
  species <- tree$tip.label
  
  
  tree <- drop.tip(tree,species[!(species %in% c("D.yakuba", "D.melanogaster", "D.simulans"))])
  tree$tip.label <- c("Dyak","Dsim","Dmel")
  a <- anc.ML(tree,vec) 
  anc <- a$ace[1]
  
  
  return(anc)
}

calc_correlations <- function(sex, species, anova, x_or_a){
  
  
  data <- get_direction_of_gexp_change(sex,species,anova,x_or_a = x_or_a)
  if(sex == 'm'){sex_long <- 'Male'}
  if(sex == 'f'){sex_long <- 'Female'}
  data <- data.frame(m = data[[1]], f = data[[2]],sex= sex_long,species = species,loc=x_or_a)
  r = cor.test(data$m,data$f)
  r_df <- data.frame(r = paste("r=", signif(r$estimate,2)),sex= sex_long,species = species,loc=x_or_a)
  p_df <- data.frame(p = paste("p=", signif(r$p.value,2)),sex= sex_long,species = species,loc=x_or_a)
  
  list(data,r_df,p_df)
}

generate_plots <- function(x_or_a){
  sexes <- c("m","f")
  species <- c("Dmel","Dsim","Dyak")
  anova <- load_anova_data()
  data <- data.frame()
  r <- data.frame()
  p <- data.frame()
  for( i in sexes){
    for(j in species){
      analysis_out <- calc_correlations(i,j,anova,x_or_a)
      data <- rbind(data,analysis_out[[1]])
      r <- rbind(r,analysis_out[[2]])
      p <- rbind(p,analysis_out[[3]])
    }
  }
  plot <- ggplot(data,aes(x=m,y=f)) + geom_point() + theme_classic()+facet_grid(species ~ sex) +
    geom_text(data = p,aes(x = 0, y = 2,label = p),hjust = 0.5,vjust = 0.5,size=3)  +
    geom_text(data = r,aes(x = 0, y = 1.5,label = r),hjust = 0.5,vjust = 0.5,size=3)+
    xlim(-1,2.7) + ylim(-1,2.7)+ xlab('Change in Male Expression') + ylab('Change in Female Expression')
  return(plot)
}

#####___Make Plots
#px <- generate_plots('X')
#pa <- generate_plots('A')

####_______________________On a gene by gene basis: Identify
####genetic constraints on sex-biased genes in Drosophila melanogaster
library(jmuOutlier)

#normalize the TPM using TMM normalization
normalize <- function(dmel_tpm){
  
  #normalize all using the TMM method from edgeR. 
  all_data <- dmel_tpm
  scale <- calcNormFactors(all_data[,(2:dim(all_data)[2])],method = "TMM")
  all_data_norm <- all_data
  for(i in c(1:length(scale))){
    all_data_norm[,(i+1)] <- all_data[,(i+1)]/scale[i]
  }
  return(all_data_norm)
}
get_strain_names <- function(data){
  out <- unique(substr(colnames(data[,c(2:dim(data)[2])]),5,8))
  return(sort(out))
}
get_med_exp <- function(data){
  strains <- get_strain_names(data)
  print(strains)
  out <- data.frame(gene_id = data$Dmel_name)
  
  colnames(data) <- c(colnames(data)[1],substr(colnames(data[,c(2:dim(data)[2])]),5,8))
  for(i in c(1:length(strains))){
    exp <- data[,colnames(data) %in% strains[i]]
    exp_median <- rowMedians(as.matrix(exp))
    exp_median <- as.data.frame(exp_median)
    colnames(exp_median) <- strains[i]
    out <- cbind(out,exp_median)
  }
  return(out)
}

cor_gene <- function(gene, data_f, data_m){
  m <- data_m[data_m$gene_id == gene, c(2:dim(data_m)[2])]
  f <- data_f[data_f$gene_id == gene, c(2:dim(data_f)[2])]
  cor <- cor(as.numeric(m), as.numeric(f))
  return(cor)
}
split_malefemale <- function(input_data){
  data <- input_data[,c(2:dim(input_data)[2])]
  sex <- substr(colnames(data),9,9)
  male_data <- cbind(data.frame(Dmel_name = input_data$genes),data[,sex == 'm'])
  female_data <- cbind(data.frame(Dmel_name = input_data$genes),data[,sex == 'f'])
  out <- list(male_data,female_data)
  names(out) <- c("male","female")
  return(out)
}




