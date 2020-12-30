library("rtracklayer")
library("GenomicFeatures")
library(plyr)
library(ggplot2)
rm(list=ls(all=TRUE)) 
setwd('/Users/skhodursky/Desktop/LabWork/Sexbias_final/Sex_Biased_Genes_DESeq/')


##define functions
get_pi_within_gene <- function(gene,pi){
  
  gene <- gen[gene]
  pi <- pi[pi %over% gene]
  data <- data.frame(start=0,pi=median(na.omit(pi$PI)))
  return(data)
}
get_pi_around_gene <- function(gene_,pi,d){
  
  gene <- gen[gene_]
  up <- flank(gene,width = d, start = T, both=F)
  down <- flank(gene, width = d, start = F, both=F)
  pi_up <- pi[pi %over% up]
  d_up <- distance(gene,pi_up)
  data_up <- data.frame(start=-d_up,pi=pi_up$PI)
  #dont include the pi within gene
  data_up <- data_up[data_up$start != 0, ]
  data_gene <- get_pi_within_gene(gene_,pi)
  pi_down <- pi[pi %over% down]
  d_down <- distance(gene,pi_down)
  data_down <- data.frame(start=d_down,pi=pi_down$PI)
  #dont include the pi within gene
  data_down <- data_down[data_down$start != 0, ]
  
  out <- rbind(data_up,data_gene,data_down)
  return(out)
}

list_to_df <- function(list){
  df <- data.frame(BIN_START= c(),pi = c())
  for(i in c(1:length(list))){
    df <- rbind(df,as.data.frame(list[[i]]))
  }
  return(df)
}

read_in_PIdata <- function(file){
  pi <- read.csv(file,header = F, stringsAsFactors = F)
  colnames(pi) <-c("chrom","start","end","PI")
  pi <- pi[(pi$start != "?") & (pi$end != "?"),]
  pi$start <- as.numeric(pi$start)
  pi$end <- as.numeric(pi$end)
  pi$PI <- as.numeric(pi$PI)
  return(makeGRangesFromDataFrame(pi,keep.extra.columns = TRUE))
}


aggregate_pi <- function(gene_list,type,chrom,pi,exp,gen_df,d = 20000){
  # PI vs distance data frame
  pi_v_d <- data.frame(start = c() , pi = c())
  for(c in chrom){
    
    if(c != '4'){
      ch <- gen_df[grepl(c,gen_df$seqnames) & gen_df$gene_id %in% exp,] 
    }else{
      ch <- gen_df[gen_df$seqnames == c & gen_df$gene_id %in% exp,] 
    }
    if(length(ch$gene_id[(ch$gene_id %in% gene_list)])>0){
      pi_v_d_chrom <- lapply(ch$gene_id[(ch$gene_id %in% gene_list)], get_pi_around_gene,pi = pi,d=d)
      pi_v_d <- rbind(pi_v_d,list_to_df(pi_v_d_chrom))
    }
  }
  pi_v_d[,"type"] <- type
  return(pi_v_d)
  
}


# function that combines pi across different classes of genes. The names of the classes should be in 
# labels_list
combine_pi_glists<- function(gene_lists, labels_list,chroms, pi, exp, gen_df){
  pi_v_d <- data.frame(start = c(), pi = c(), type = c())  
  for(i in seq(length(gene_lists))){
    pi_v_d <- rbind(pi_v_d,
                    aggregate_pi(gene_lists[[i]], labels_list[[i]], chroms, pi, exp, gen_df))
  }
  return(pi_v_d)
}
