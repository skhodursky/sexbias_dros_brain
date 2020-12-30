library(GenomicFeatures)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
setwd('/Users/skhodursky/Desktop/LabWork/Sexbias_final/Analysis/')

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

adjust_pval <- function(pval_vector){
  p <- p.adjust(c(pval_vector),method = "BH" )
  p <- as.character(symnum(p, corr = FALSE, na = FALSE,
                           cutpoints = c(0,1e-5, 1e-4, 1e-3, 1e-2, 5e-2,1),
                           symbols = c("*****","****","***","**","*",'NS')))
  return(p)
}



#function to test for the significance of MS/MSR using simulation
# stat is the statistic of interest (i.e. MS_sp, MS_st etc.)
# geneset is geneset of interest (i.e. male-biased X chromsome genes)
get_pvalue_sim <- function(anova_data, geneset, stat, numsim){
  #define sampling function
  sample_genomewide <- function(iter,anova_data, n){
    idx <- sample(c(1:dim(anova_data)[1]),n)
    return(median(anova_data[idx,stat]))
  }
  #get number of genes sampled and the observed median F-ratio
  num_sampled <- sum(geneset %in% anova_data$gene_id)
  median_ratio_sampled <- median(anova_data[anova_data$gene_id %in% geneset, stat])
  
  iteration <- c(1:numsim)
  sim_out <- sapply(iteration, sample_genomewide, anova_data = anova_data, n = num_sampled)
  
  if(median_ratio_sampled > median(anova_data[,stat])){
    p <- 2*(sum(sim_out > median_ratio_sampled)/numsim)
  }else{
    p <- 2*(sum(sim_out < median_ratio_sampled)/numsim)
  }
  return(p)
}

# the input anova_data should be log transformed.
simsig_plots <- function(stat, anova_mal, anova_fem, numsim = 1000000){

  
  #Load gene lists
  sex_biased_lists <- load_sbg_data()
  x_m <-sex_biased_lists$x_m
  x_f <-sex_biased_lists$x_f
  auto_m <-sex_biased_lists$auto_m
  auto_f <-sex_biased_lists$auto_f
  x_u <-sex_biased_lists$x_u
  auto_u <-sex_biased_lists$auto_u
  
  #get simulation p-values for each geneset
  genelists <- list(x_f, x_m, x_u,
                    auto_f, auto_m, auto_u)
  p_mal <- c()
  p_fem <- c()
  for(geneset in genelists){
    p_mal <- c(p_mal,get_pvalue_sim(anova_mal, geneset, stat, numsim))
    p_fem <- c(p_fem,get_pvalue_sim(anova_fem, geneset, stat, numsim))
    
  }
  
  
  
  #multiple test adjust p-values
  pvals <- data.frame(p = adjust_pval(c(p_mal, p_fem)))
  #make data_frames for plotting
  anova_m <- rbind(data.frame(Stat = anova_mal[anova_mal$gene_id %in% x_f,stat], class = 'Female X' ),
                   data.frame(Stat = anova_mal[anova_mal$gene_id %in% x_m,stat], class = 'Male X' ),
                   data.frame(Stat= anova_mal[anova_mal$gene_id %in% x_u,stat], class = 'Unb. X'),
                   data.frame(Stat = anova_mal[anova_mal$gene_id %in% auto_f,stat], class = 'Female Auto.' ),
                   data.frame(Stat= anova_mal[anova_mal$gene_id %in% auto_m,stat], class = 'Male Auto.'),
                   data.frame(Stat= anova_mal[anova_mal$gene_id %in% auto_u,stat], class = 'Unb. Auto.'))
  
  anova_f <- rbind(data.frame(Stat = anova_fem[anova_fem$gene_id %in% x_f,stat], class = 'Female X'),
                   data.frame(Stat = anova_fem[anova_fem$gene_id %in% x_m,stat], class = 'Male X'),
                   data.frame(Stat= anova_fem[anova_fem$gene_id %in% x_u,stat], class = 'Unb. X'),
                   data.frame(Stat= anova_fem[anova_fem$gene_id %in% auto_f,stat], class = 'Female Auto.'),
                   data.frame(Stat= anova_fem[anova_fem$gene_id %in% auto_m,stat], class = 'Male Auto.'),
                   data.frame(Stat= anova_fem[anova_fem$gene_id %in% auto_u,stat], class = 'Unb. Auto.'))
  
  
  anova_m$class <- factor(anova_m$class)
  anova_f$class <- factor(anova_f$class)
  
  
  #Plot details
  if(stat == 'MSR'){
    ylabel <- expression('Log('*italic('MS'['sp']/'MS'['st']) * ')')}
  if(stat == 'MS_sp'){
    ylabel <- expression('Log('*italic('MS'['sp']) * ')')}
  if(stat == 'MS_st'){
    ylabel <- expression('Log('*italic('MS'['st']) * ')')}
  plot_bot <- min(c(anova_m$Stat, anova_f$Stat))
  plot_top <-  max(c(anova_m$Stat, anova_f$Stat))
  
  
  y_position <- plot_top + 1
  
  ##plot for genes in males
  plot_m <- ggplot(anova_m,aes(x=class,y=(Stat),fill = class)) + geom_boxplot()  + 
    annotate(geom="text",x=c(1:6),y=y_position,label=pvals$p[c(1:6)], size = 4) +
    geom_hline(yintercept = median(anova_m$Stat), linetype="dashed") + 
    xlab('') + ylim(plot_bot,plot_top + 2) + ylab(ylabel)+
    ggtitle('Expression in males') +
    theme_classic()+ theme(legend.position="none")+
    scale_fill_brewer(palette="Dark2")
  
  
  
  plot_f <- ggplot(anova_f,aes(x=class,y=(Stat),fill = class)) + geom_boxplot()  + 
    annotate(geom="text",x=c(1:6),y=y_position,label=pvals$p[c(7:12)], size = 4) +
    geom_hline(yintercept = median(anova_f$Stat), linetype="dashed") + 
    xlab('') + ylim(plot_bot,plot_top + 2) + ylab(ylabel)+
    ggtitle('Expression in females') +
    theme_classic()+ theme(legend.position="none")+
    scale_fill_brewer(palette="Dark2")
  
  plot <- plot_grid(plot_f+theme(text = element_text(size=10.5),plot.title=element_text(size=10))+ 
                      scale_x_discrete(labels=c("Female\nX","Male\nX", 
                                                "Unbiased\nX", "Female\nA",
                                                "Male\nA","Unbiased\nA")),
                    plot_m+theme(text = element_text(size=10.5),plot.title=element_text(size=10))+
                      scale_x_discrete(labels=c("Female\nX","Male\nX", "Unbiased\nX",
                                                "Female\nA","Male\nA","Unbiased\nA")),
                    nrow = 1,align = T)
  return(plot)
}

#this is to test for enrichment of stabilizing or relaxed selection in genes. 
test_stabilizing <- function(anova_data, geneset){
    
  ms_sp_sig <- anova_data$gene_id[anova_data$MS_sp < median(anova_data$MS_sp)]
  ms_st_sig <- anova_data$gene_id[anova_data$MS_st < median(anova_data$MS_st)]
  #total number of genes with MS_sp and MS_st less than genome wide medians
  num_sig <- sum(ms_sp_sig %in% ms_st_sig )
  #number of genes in geneset found in anova dataframe
  num_sampled <- sum(geneset %in% anova_data$gene_id)
  #number of genews with MS_sp and MS_st less than genome wide medians in geneset
  num_sig_geneset <- sum((geneset %in% ms_sp_sig) & (geneset %in% ms_st_sig) )
  #number of genes in anova dataframe
  total <- dim(anova_data)[1]
  #fold enrichment
  fe <- (num_sig_geneset/num_sampled)/(num_sig/total)
  #multiply the p value by two because we're changing the tail for the test
  if(fe>=1){
    p <- 2*phyper(num_sig_geneset, num_sig, total-num_sig,num_sampled, lower.tail = FALSE)
  }else{
    p <- 2*phyper(num_sig_geneset, num_sig, total-num_sig,num_sampled, lower.tail = TRUE)
  }
  out <- list(p, fe)
  names(out) <- c('p','fe')
  return(out)
}

stabilizing_plots <- function(anova_mal, anova_fem){
  library(ggplot2)
  library(cowplot)
  library(RColorBrewer)
  
  
  
  #Load gene lists
  sex_biased_lists <- load_sbg_data()
  x_m <-sex_biased_lists$x_m
  x_f <-sex_biased_lists$x_f
  auto_m <-sex_biased_lists$auto_m
  auto_f <-sex_biased_lists$auto_f
  x_u <-sex_biased_lists$x_u
  auto_u <-sex_biased_lists$auto_u
  
  #get simulation p-values for each geneset
  genelists <- list(x_f, x_m, x_u,
                    auto_f, auto_m, auto_u)
  p_mal <- c()
  p_fem <- c()
  for(geneset in genelists){
    p_mal <- c(p_mal,test_stabilizing(anova_mal, geneset)$p)
    p_fem <- c(p_fem,test_stabilizing(anova_fem, geneset)$p)
    
  }
  
  
  
  #multiple test adjust p-values
  pvals <- data.frame(p = adjust_pval(c(p_mal,p_fem)))
  #make data_frames for plotting
  anova_m <- rbind(data.frame(FE = test_stabilizing(anova_mal, x_f)$fe, class = 'Female X' ),
                   data.frame(FE = test_stabilizing(anova_mal, x_m)$fe, class = 'Male X' ),
                   data.frame(FE= test_stabilizing(anova_mal, x_u)$fe, class = 'Unb. X'),
                   data.frame(FE = test_stabilizing(anova_mal, auto_f)$fe, class = 'Female Auto.' ),
                   data.frame(FE= test_stabilizing(anova_mal, auto_m)$fe, class = 'Male Auto.'),
                   data.frame(FE= test_stabilizing(anova_mal, auto_u)$fe, class = 'Unb. Auto.'))
  
  anova_f <- rbind(data.frame(FE = test_stabilizing(anova_fem, x_f)$fe, class = 'Female X' ),
                   data.frame(FE = test_stabilizing(anova_fem, x_m)$fe, class = 'Male X' ),
                   data.frame(FE= test_stabilizing(anova_fem, x_u)$fe, class = 'Unb. X'),
                   data.frame(FE = test_stabilizing(anova_fem, auto_f)$fe, class = 'Female Auto.' ),
                   data.frame(FE= test_stabilizing(anova_fem, auto_m)$fe, class = 'Male Auto.'),
                   data.frame(FE= test_stabilizing(anova_fem, auto_u)$fe, class = 'Unb. Auto.'))
  
  
  anova_m$class <- factor(anova_m$class)
  anova_f$class <- factor(anova_f$class)
  
  
  #Plot details
  
  ylabel <- expression(atop('FE of genes with both low' 
                         , italic('MS'['sp']) ~ 'and low' ~ italic('MS'['st'])))
  plot_bot <- 0
  plot_top <-  max(c(anova_m$FE, anova_f$FE))
  
  
  y_position <- plot_top + .1
  
  ##plot for genes in males
  plot_m <- ggplot(anova_m, aes(x=class, y=FE, fill = class)) + geom_col()  + 
    annotate(geom="text", x=c(1:6), y=y_position, label=pvals$p[c(1:6)], size = 4) +
    geom_hline(yintercept = 1,  linetype="dashed") + 
    xlab('') + ylim(plot_bot, plot_top + .2) + ylab(ylabel)+
    ggtitle('Expression in males') +
    theme_classic()+ theme(legend.position="none")+
    scale_fill_brewer(palette="Dark2")
  
  
  # plot for genes in females
  plot_f <- ggplot(anova_f, aes(x=class, y=FE, fill = class)) + geom_col()  + 
    annotate(geom="text", x=c(1:6), y=y_position, label=pvals$p[c(7:12)], size = 4) +
    geom_hline(yintercept = 1,  linetype="dashed") + 
    xlab('') + ylim(plot_bot, plot_top + .2) + ylab(ylabel)+
    ggtitle('Expression in females') +
    theme_classic()+ theme(legend.position="none")+
    scale_fill_brewer(palette="Dark2")
  
  plot <- plot_grid(plot_f+theme(text = element_text(size=10.5), plot.title=element_text(size=10))+ 
                      scale_x_discrete(labels=c("Female\nX", "Male\nX",  "Unbiased\nX",  "Female\nA", "Male\nA", "Unbiased\nA")), 
                    plot_m+theme(text = element_text(size=10.5), plot.title=element_text(size=10))+
                      scale_x_discrete(labels=c("Female\nX", "Male\nX",  "Unbiased\nX",  "Female\nA", "Male\nA", "Unbiased\nA")), 
                    nrow = 1, align = T)
  return(plot)
  
}


simsig_newold_plots <- function(stat, anova_mal, anova_fem, 
                                new_m, new_f, old_m, old_f, numsim = 1000000){
  library(ggplot2)
  library(cowplot)
  library(RColorBrewer)
  
  
  
  #Load and define gene lists
  sex_biased_lists <- load_sbg_data()
  x_m <-sex_biased_lists$x_m
  x_f <-sex_biased_lists$x_f
  auto_m <-sex_biased_lists$auto_m
  auto_f <-sex_biased_lists$auto_f
 
  x_m_n <-x_m[x_m %in% new_m]
  x_m_o <-x_m[x_m %in% old_m]
  
  x_f_n <-x_f[x_f %in% new_f]
  x_f_o <-x_f[x_f %in% old_f]
   
  auto_f_n <-auto_f[auto_f %in% new_f]
  auto_f_o <-auto_f[auto_f %in% old_f]
  
  auto_m_n <-auto_m[auto_m %in% new_m]
  auto_m_o <-auto_m[auto_m %in% old_m]
  
  #get simulation p-values for each geneset
  genelists <- list(x_f_n,x_f_o,x_m_n,x_m_o,
                    auto_f_n,auto_f_o,auto_m_n,auto_m_o)
  p_mal <- c()
  p_fem <- c()
  for(geneset in genelists){
    p_mal <- c(p_mal,get_pvalue_sim(anova_mal, geneset, stat,numsim))
    p_fem <- c(p_fem,get_pvalue_sim(anova_fem, geneset, stat,numsim))
    
  }
  
  
  
  #multiple test adjust p-values
  pvals <- data.frame(p = adjust_pval(c(p_mal,p_fem)))
  #make data_frames for plotting
  anova_m <- rbind(data.frame(Stat = anova_mal[anova_mal$gene_id %in% x_f_n,stat], class = 'xfn' ),
                   data.frame(Stat = anova_mal[anova_mal$gene_id %in% x_f_o,stat], class = 'xfo' ),
                   data.frame(Stat = anova_mal[anova_mal$gene_id %in% x_m_n,stat], class = 'xmn'),
                   data.frame(Stat = anova_mal[anova_mal$gene_id %in% x_m_o,stat], class = 'xmo' ),
                   data.frame(Stat = anova_mal[anova_mal$gene_id %in% auto_f_n,stat], class = 'afn'),
                   data.frame(Stat = anova_mal[anova_mal$gene_id %in% auto_f_o,stat], class = 'afo'),
                   data.frame(Stat = anova_mal[anova_mal$gene_id %in% auto_m_n,stat], class = 'amn'),
                   data.frame(Stat = anova_mal[anova_mal$gene_id %in% auto_m_o,stat], class = 'amo'))
  
  anova_f <- rbind(data.frame(Stat = anova_fem[anova_fem$gene_id %in% x_f_n,stat], class = 'xfn' ),
                   data.frame(Stat = anova_fem[anova_fem$gene_id %in% x_f_o,stat], class = 'xfo' ),
                   data.frame(Stat = anova_fem[anova_fem$gene_id %in% x_m_n,stat], class = 'xmn'),
                   data.frame(Stat = anova_fem[anova_fem$gene_id %in% x_m_o,stat], class = 'xmo' ),
                   data.frame(Stat = anova_fem[anova_fem$gene_id %in% auto_f_n,stat], class = 'afn'),
                   data.frame(Stat = anova_fem[anova_fem$gene_id %in% auto_f_o,stat], class = 'afo'),
                   data.frame(Stat = anova_fem[anova_fem$gene_id %in% auto_m_n,stat], class = 'amn'),
                   data.frame(Stat = anova_fem[anova_fem$gene_id %in% auto_m_o,stat], class = 'amo'))
  
  
  anova_m$class <- factor(anova_m$class)
  anova_f$class <- factor(anova_f$class)
  
  
  #Plot details
  if(stat == 'MSR'){
    ylabel <- expression('Log('*italic('MS'['sp']/'MS'['st']) * ')')}
  if(stat == 'MS_sp'){
    ylabel <- expression('Log('*italic('MS'['sp']) * ')')}
  if(stat == 'MS_st'){
    ylabel <- expression('Log('*italic('MS'['st']) * ')')}
  plot_bot <- min(c(anova_m$Stat, anova_f$Stat))
  plot_top <-  max(c(anova_m$Stat, anova_f$Stat))
  
  
  y_position <- plot_top + 1
  
  ##plot for genes in males
  plot_m <- ggplot(anova_m,aes(x=class,y=(Stat),fill = class)) + geom_boxplot()  + 
    annotate(geom="text",x=c(1:8),y=y_position,label=pvals$p[c(1:8)], size = 4) +
    geom_hline(yintercept = median(anova_m$Stat), linetype="dashed") + 
    xlab('') + ylim(plot_bot,plot_top + 2) + ylab(ylabel)+
    ggtitle('Expression in males') +
    theme_classic()+ theme(legend.position="none")+
    scale_fill_brewer(palette="Dark2")
  
  
  
  plot_f <- ggplot(anova_f,aes(x=class,y=(Stat),fill = class)) + geom_boxplot()  + 
    annotate(geom="text",x=c(1:8),y=y_position,label=pvals$p[c(9:16)], size = 4) +
    geom_hline(yintercept = median(anova_f$Stat), linetype="dashed") + 
    xlab('') + ylim(plot_bot,plot_top + 2) + ylab(ylabel)+
    ggtitle('Expression in females') +
    theme_classic()+ theme(legend.position="none")+
    scale_fill_brewer(palette="Dark2")
  
  plot <- plot_grid(plot_f+theme(text = element_text(size=10.5),plot.title=element_text(size=10))+ 
                      scale_x_discrete(labels=c("F, 1sp.\nX","F, 3sp.\nX" ,"M, 1sp.\nX","M, 3sp.\nX" ,
                                                "F, 1sp.\nA","F, 3sp.\nA", "M, 1sp.\nA","M, 3sp.\nA" )),
                    plot_m+theme(text = element_text(size=10.5),plot.title=element_text(size=10))+
                      scale_x_discrete(labels=c("F, 1sp.\nX","F, 3sp.\nX" ,"M, 1sp.\nX","M, 3sp.\nX" ,
                                                "F, 1sp.\nA","F, 3sp.\nA" , "M, 1sp.\nA","M, 3sp.\nA" )),
                    nrow = 1,align = T)
  return(plot)
}
