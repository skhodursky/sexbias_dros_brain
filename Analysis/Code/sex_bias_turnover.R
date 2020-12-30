rm(list=ls(all=TRUE))
library(GenomicFeatures)


# get the top n sex-biased genes in each species
get_topN_sbgenes <- function(species, num){
  deseq_summary <- read.csv(paste0(species, '/total_DESeq_output.txt'))
  deseq_summary <- deseq_summary[order(deseq_summary$pvalue),]
  below_cutoff <- deseq_summary[c(1:num),]
  return(below_cutoff)
}

# add dmel orthologs to dyak and dsim
add_dmel_ortho <- function(species, data, orthologs){
  dmel_ortho <- function(gene, species, orthologs){
    return(orthologs[orthologs[, species] == gene, 'Dmel'])
  }
  mel_gene <- sapply(data$gene_id, dmel_ortho, species = species, orthologs = orthologs)
  data <- cbind(data.frame(Dmel_name= mel_gene), data)
  return(data)
}



# get ovelaps between sets
get_overlaps <- function(a, b, c){
  # get genes in all 3, 2, and only 1 species
  in_3 <- Reduce(intersect, list(a, b, c))
  in_2 <- unique(c(setdiff(intersect(a, b), c),
                   setdiff(intersect(a, c), b), 
                   setdiff(intersect(b, c), a)))
  in_1 <- unique(c(setdiff(c(a, b, c), c(in_2, in_3))))
  
  return(list(in_1,in_2,in_3))
  
}

analyze_turnover <- function(){
  orthologs <- read.csv('Evolution_of_Sex_Bias/one_one_orthologs_mel_yak_sim.csv', stringsAsFactors = FALSE)
  
  sbg_mel <- read.csv('Dmel/sex_biased_genes_dmel.txt', stringsAsFactors = FALSE)
  #get number of sex-biased genes in melanogaster
  num <- dim(sbg_mel)[1]
  sbg_mel <- sbg_mel[,c("gene_id", "sex")]
  sbg_mel <- sbg_mel[sbg_mel$gene_id %in% orthologs$Dmel,]
  colnames(sbg_mel) <- c("Dmel_name","sex")
  male_mel <- as.vector(sbg_mel$Dmel_name[sbg_mel$sex == 'm'])
  female_mel <- as.vector(sbg_mel$Dmel_name[sbg_mel$sex == 'f'])
  
  # get top genes in dsim
  sbg_sim <- get_topN_sbgenes('dsim', num)
  sbg_sim <- sbg_sim[,c("gene_id", "sex")]
  sbg_sim <- sbg_sim[sbg_sim$gene_id %in% orthologs$Dsim,]
  sbg_sim <- add_dmel_ortho('Dsim', sbg_sim, orthologs)
  male_sim <- as.vector(sbg_sim$Dmel_name[sbg_sim$sex == 'm'])
  female_sim <- as.vector(sbg_sim$Dmel_name[sbg_sim$sex == 'f'])
  
  # get top genes in dyak
  sbg_yak <- get_topN_sbgenes('dyak', num)
  sbg_yak <- sbg_yak[,c("gene_id", "sex")]
  sbg_yak <- sbg_yak[sbg_yak$gene_id %in% orthologs$Dyak,]
  sbg_yak <- add_dmel_ortho('Dyak', sbg_yak, orthologs)
  male_yak <- as.vector(sbg_yak$Dmel_name[sbg_yak$sex == 'm'])
  female_yak <- as.vector(sbg_yak$Dmel_name[sbg_yak$sex == 'f'])
  
  #remove switched genes before looking at species overlaps
  male_all <- unique(c(male_mel,male_sim,male_yak))
  female_all <- unique(c(female_mel,female_sim,female_yak))
  ##____________look at switched genes
  switched <- male_all[male_all %in% female_all]
  switched_df <- data.frame(gene_id = switched,dmel_sex ='NS',
                            dsim_sex= 'NS', dyak_sex= 'NS', 
                            stringsAsFactors = FALSE)
  switched_df$dmel_sex[switched_df$gene_id %in% male_mel] <- 'm'
  switched_df$dmel_sex[switched_df$gene_id %in% female_mel] <- 'f'
  
  switched_df$dsim_sex[switched_df$gene_id %in% male_sim] <- 'm'
  switched_df$dsim_sex[switched_df$gene_id %in% female_sim] <- 'f'
  
  switched_df$dyak_sex[switched_df$gene_id %in% male_yak] <- 'm'
  switched_df$dyak_sex[switched_df$gene_id %in% female_yak] <- 'f'
  # seperate male and female genes and find set overlaps
  male_mel <- setdiff(male_mel, switched)
  female_mel <- setdiff(female_mel, switched)
  
  male_sim <- setdiff(male_sim, switched)
  female_sim <- setdiff(female_sim, switched)
  
  male_yak <- setdiff(male_yak, switched)
  female_yak <- setdiff(female_yak, switched)
  
  #look at species overlaps for male and female-biased genes
  sets_male <- get_overlaps(male_mel, male_sim, male_yak)
  sets_female <- get_overlaps(female_mel, female_sim, female_yak)
  
  out <- list('male_1' = sets_male[[1]],'male_2' = sets_male[[2]],'male_3' = sets_male[[3]],
              'female_1' = sets_female[[1]],'female_2' = sets_female[[2]],'female_3' = sets_female[[3]],'switched' = switched )
  return(out)
}
