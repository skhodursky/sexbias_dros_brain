library("rtracklayer")
library("GenomicFeatures")
library("IRanges")
library(plyr)
rm(list=ls(all=TRUE)) 
#setwd('/Users/skhodursky/Desktop/LabWork/Sexbias_final/Sex_Biased_Genes_DESeq/An')


##define functions 

#This function gets you the nucleotide diversity (hk) around a given gene d bases upstream and downstream
#from gene not including the gene itself
get_hk_around_gene <- function(gene,hk,d = 20000,gen){
  g <- gen[gene,]
  #exclude region containing gene
  hk <- hk[!(hk %over% g),]
  g_up <- flank(g,d,start = T, both = F)
  g_down <- flank(g,d,start = F, both = F)
  score_up <- min(na.omit(hk[hk %over% g_up,]$score))
  score_down <- min(na.omit(hk[hk %over% g_down,]$score))
  
  minimum <- min(score_up,score_down)
  if(is.finite(minimum)){
    return(minimum)
  }else{
    return(NA)
  }
  
}


##loads hk data. Make sure data is in appropriate coordinate system (version 6)
read_in_hkdata <- function(file){
  hk <- import.bedGraph(file)
  return(hk)
}

p_to_star <- function(p){
  p <- as.character(symnum(p, corr = FALSE, na = FALSE,
                           cutpoints = c(0,1e-5, 1e-4, 1e-3, 1e-2, 5e-2,1),
                           symbols = c("*****","****","***","**","*",'NS')))
  return(p)
}

# aggregate hkal data

combine_hkal_glists<- function(gene_lists, labels_list, hk, gen){
  hkal_df<- data.frame(hkal = c(), type = c())  
  for(i in seq(length(gene_lists))){
    hkal_df <- rbind(hkal_df,
                     data.frame(hkal = sapply(gene_lists[[i]],get_hk_around_gene,hk=hk,gen=gen ),
                                type=labels_list[[i]]))
  }
  return(hkal_df)
}
