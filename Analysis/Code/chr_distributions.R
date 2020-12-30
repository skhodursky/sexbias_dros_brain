library(ggplot2)
library(gridExtra)
library(cowplot)
rm(list=ls(all=TRUE)) 

# load sex-biased genes
load_sbdata <- function(species){
  sbg_df <- read.csv(paste0(species, '/sex_biased_genes_', species, '.txt'), stringsAsFactors = FALSE)
  ref <- read.csv(paste0(species, '/genes_with_finite_p_val.txt'), stringsAsFactors = FALSE)
  # just to be safe
  ref <- ref[!duplicated(ref$gene_id),]
  out <- list(sbg_df,ref)
  names(out) <- c("sbg_df","ref")
  return(out)
}

# randomly sample genes from pool of genes with finite p-values according to deseq
simulate <- function(i, sex, sbg_df, ref){
  # number of observed sb genes
  num <- sum(sbg_df$sex == sex)
  # random sample of "expressed" genes
  s <- sample(c(1:dim(ref)[1]), num, replace = F)
  h <- ref[s,]$seqnames
  # count the number of sampled genes on each chromosome
  nX <- sum(grepl('X', h))
  n2L <- sum(grepl('2L', h))
  n2R <- sum(grepl('2R', h))
  n3L <- sum(grepl('3L', h))
  n3R <- sum(grepl('3R', h))
  distribution <- data.frame(nX = nX, n2L = n2L, n2R= n2R, n3L= n3L, n3R= n3R)
  return(distribution)
}

# this is effectively a 2 tailed test for observed relative to simulated distributions
get_pval <- function(chr, tot_distribution, sbg_df, sex){
  # count the observed number of sex-biased genes on given chr
  #num <- sum((sbg_df$sex == sex) & grepl(chr, sbg_df$seqnames))
  num <- get_numobs(chr, sbg_df, sex)
  chr <- paste0('n', chr)
  
  if(num < median(tot_distribution[, chr])){
    out <- 2*sum(tot_distribution[, chr] <= num)/length(tot_distribution[,chr] )
  }
  else if(num >= median(tot_distribution[, chr])){
    out <- 2*sum(tot_distribution[, chr] >= num)/length(tot_distribution[,chr] )
  }  
  out <- as.character(symnum(c(out), corr = FALSE, na = FALSE,
                             cutpoints = c(0,1e-5, 1e-4, 1e-3, 1e-2, 5e-2,1),
                             symbols = c("*****","****","***","**","*",'')))
  
  return(out)
}

# count the number of sex-biased genes on a chromosome
get_numobs <- function(chr, sbg_df, sex){
  # count the observed number of sex-biased genes on given chr
  num <- sum((sbg_df$sex == sex) & grepl(chr, sbg_df$seqnames))
  return(num)
}

# function to reshape dataframe
reshape_dist_df <- function(tot_distribution){
  out <- data.frame()
  cols <- colnames(tot_distribution)
  for(c in cols){
    number <- data.frame(counts = tot_distribution[, c],chrom = substr(c, 2, nchar(c)))
    out <- rbind(out, number)
  }
  return(out)
}
# plot simulation results
plot_simulations <- function(numsim, sex, species){
  
  sbg_df <- load_sbdata(species)$sbg_df
  ref <- load_sbdata(species)$ref
  
  tot_distribution <- do.call(rbind,lapply((c(1:numsim)), simulate,
                                           sbg_df=sbg_df, sex=sex, ref=ref))
  
  # for each chromosome get the number of sex-biased genes observed
  # and p-value for that observation
  chroms <- c('X', '2L', '2R', '3L', '3R')
  obs_num <- sapply(chroms, get_numobs, sbg_df = sbg_df, sex = sex)
  p_vals <- sapply(chroms, get_pval, tot_distribution = tot_distribution, 
                   sbg_df = sbg_df, sex = sex)
  real_counts <- data.frame(chrom = chroms, rc = obs_num)
  pvals <- data.frame(chrom = chroms, p = p_vals)
  
  # specify y labels
  if (sex == 'm'){yl <- 'Male'}
  if(sex == 'f'){yl <- 'Female'}
  ylimit <- max(obs_num)*1.2
  plot <- ggplot(reshape_dist_df(tot_distribution), aes(y= counts)) + theme_classic() + 
    geom_boxplot(width = 0.1) + geom_point(data=real_counts, aes(x=0, y=rc), size=2, color = "red") + 
    ylab('Number of Genes') +ggtitle(paste(yl, '-biased genes', sep='')) + ylim(0, ylimit) + 
    geom_text(data = pvals, aes(x = 0, y = 10, label = p), hjust = 0.5, vjust = 0.5, size=4)  +
    facet_grid(. ~ chrom) + xlab('') + 
    theme( legend.position = "none", axis.text.x=element_blank(), 
           axis.ticks.x=element_blank(), plot.title = element_text(hjust = 0.5))
  
  return(plot)
}
