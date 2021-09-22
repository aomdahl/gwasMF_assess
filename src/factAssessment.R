####### Assess the quality of a factorization
####### Made by Ashton Omdahl, Battle Lab, Aug 2021
rm(list = ls()) #clear the environment
pacman::p_load(data.table, tidyr, dplyr, readr, ggplot2, stringr, Xmisc, PRROC, ROCR,svMisc, stats, cowplot)
source("/work-zfs/abattle4/ashton/snp_networks/scratch/ldsc_all_traits/src/assessment_functions.R")
#######functions


################### MAIN ###################

parser <- ArgumentParser$new()
parser$add_description("Script to assess how well a matrix factorization captures tissue effects. Performs in both trait and factor specific manner.")
parser$add_argument("--factors", type = 'character', help = "What is the factorization file" )
parser$add_argument("--output", type = 'character', help = "output directory to write to")
parser$add_argument("--ldsc_reference", type = 'character', help = "This directory contains the true LDSC information for each trait")
parser$add_argument("--ldsc_dir", type = 'character', default = "/work-zfs/abattle4/ashton/snp_networks/gwas_decomp_ldsc/results/ldsc/seed1_thresh0.7/", help = "Where is the ldsc output for each factor")
parser$add_argument("--trait.ids", type = 'character', default = "/work-zfs/abattle4/ashton/snp_networks/gwas_decomp_ldsc/trait_selections/seed2_thresh0.9_h2-0.1.studies.tsv", help = "Where are the ids for the traits?")
parser$add_argument("--trait.names", type = 'character', default = "/work-zfs/abattle4/ashton/snp_networks/gwas_decomp_ldsc/trait_selections/seed2_thresh0.9_h2-0.1.names.tsv", help = "Where are the names for the traits?")
parser$add_argument("--n_overlap", type = 'integer', default = 2, help = "overlap in tissues.")
parser$add_argument("--all", type = 'logical', default = FALSE, action='store_true', help = "Run all analysis metrics.")
parser$add_argument("--trait_specific", type = 'logical', default = FALSE, action='store_true', help = "Run trait_specific analysis metrics.")
parser$add_argument("--factor_specific", type = 'logical', default = FALSE, action='store_true', help = "Run factor-specific analysis metrics.")
parser$add_argument("--simple", type = 'logical', default = FALSE, help = "Run simple tests.",action='store_true')
parser$add_argument("--fdr.thresh", type = 'double', default = 0.05, help = "True value cutoff.")
parser$add_argument("--overlap_test", type = 'logical', action='store_true', default = FALSE, help = "Run overlap test.")
parser$add_argument('--help',type='logical',action='store_true',help='Print the help page', default = FALSE)
parser$helpme()
args <- parser$get_args()
if(FALSE)
{
  T="/work-zfs/abattle4/ashton/snp_networks/custom_l1_factorization/"
  args <- list()
  #args$factors <- paste0("/work-zfs/abattle4/ashton/snp_networks/custom_l1_factorization/factorization_data/beta_se_fixed_first900_900.factors.txt")
  #args$factors <- paste0("/work-zfs/abattle4/ashton/snp_networks/custom_l1_factorization/factorization_data/aug_28.ashr.factors.txt")
  args$factors <- paste0("/work-zfs/abattle4/ashton/snp_networks/custom_l1_factorization/factorization_data/A400_L400_B_SE.factors.txt")
  plotFactors(fread(args$factors), trait.names, title = "A400_L400")
  args$ldsc_reference <-"/work-zfs/abattle4/ashton/snp_networks/scratch/ldsc_all_traits/ldsc_results/seed2_thres0.9_h2-0.1/"
  #args$ldsc_dir <-paste0(T,"/results/beta_se_fixed_first900_900/ldsc_enrichment_Multi_tissue_chromatin/" ) 
  args$ldsc_dir <-"/work-zfs/abattle4/ashton/snp_networks/custom_l1_factorization/results/aug_28.ashr//ldsc_enrichment_Multi_tissue_chromatin/"
  
  args$trait.ids <- "/work-zfs/abattle4/ashton/snp_networks/gwas_decomp_ldsc/trait_selections/seed2_thresh0.9_h2-0.1.studies.tsv"
  args$trait.names <- "/work-zfs/abattle4/ashton/snp_networks/gwas_decomp_ldsc/trait_selections/seed2_thresh0.9_h2-0.1.names.tsv"
  args$n_overlap <- 2
  args$fdr.thresh <- 0.05
}


#Read in
ldsc.reference <- readInLDSC(args$ldsc_reference)

factorization.enrichment <- readInLDSC(args$ldsc_dir, type = "factorization enrichment")

factorization <- fread(args$factors)
sparsity.threshold <- max(abs(factorization)) * 1e-3 #kind of arbitrary, but at least scaled to the scale of the matrix.
trait.names <- scan("/work-zfs/abattle4/ashton/snp_networks/gwas_decomp_ldsc/trait_selections/seed2_thresh0.9_h2-0.1.names.tsv", what = character())
trait.studies <- scan("/work-zfs/abattle4/ashton/snp_networks/gwas_decomp_ldsc/trait_selections/seed2_thresh0.9_h2-0.1.studies.tsv", what = character())
K <- ncol(factorization)
ntraits <- nrow(factorization)
res <- matrix(NA, ntraits, 2)
print("Information loaded")
################### TRAIT-specific metrics ###################
if(args$all || args$trait_specific)
{
  #assess the P/R curve for each trait
  print("Creating the PR curves now....")
  clean.names <- gsub(pattern = ",", replacement = "", x = trait.names)
  for(i in 1:ntraits)
  {
    progress(i, progress.bar = TRUE)
    trait.id <- trait.studies[i]
    pr <- areaPRPlot(i, trait.studies, factorization, factorization.enrichment, ldsc.reference)
    res[i,1] <- pr$auc.integral
    res[i,2] <- pr$auc.davis.goadrich
    #save the plot:
    if(pr$auc.integral != 0)
    {
      png(paste0(args$output, "/pr_", clean.names[i], "_", trait.id, ".png") )
      plot(pr)
      dev.off() 
    } else{
      print(paste0("Unable to print for trait ", trait.names[i], " had a 0 AUPRC"))
    }
  }
  print("Successfully wrote out all the trait-specific curves.")
  
  #zero out the NAs:
  res[is.na(res[,1]), 1] <- 0
  #save a tabular version for each trait, and an image too.
  #AUPR per trait:
  print("Reporting the AUPR per trait...")
  by.trait <- data.frame("trait" = as.character(trait.names), "id" = as.character(trait.studies), 'aupr' = res[,1]) %>%
    add_row(trait = "Total_avg", id = "NA", 'aupr' = mean(res[,1]))
  write_tsv(x= by.trait, file = paste0(args$output, "/per_trait_aupr.tsv"))
  
  #image version
  plot <- data.frame("names" = trait.names,"auprc" = res[,1])
  p <- ggplot(dat = plot, aes(x = reorder(names, -auprc), y = auprc)) + geom_bar(stat = "identity")  + theme_minimal(15) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  + scale_x_discrete(label = function(x) abbreviate(x, minlength = 20)) +
    xlab("Trait") + ylab("AUPR") + geom_hline(yintercept = mean(plot$auprc), color="red")
  ggsave(filename =paste0(args$output, "/pr_aupr_traits.png"), plot = p, width = 10)
  #write_tsv(x = plot, file = paste0(args$output, "/pr_aupr_traits.txt"))
  
}

################### FACTOR-specific metrics ###################
#AUPR per factor
if(args$all || args$factor_specific)
{
  #AUPR per factor by sharing
  #####HERE 
  #Number of shared enrichments within a factor..... (how good is the grouping?)
  l <- list()
  for(i in 1:ncol(factorization))
  {
    l[[i]] <- factorSpecificAUPRShared(i, trait.studies, factorization, factorization.enrichment, ldsc.reference, sparsity.threshold, answer_thres = 0.1)
    print(i)
  } 
  num.factors.with.enrichemnts <- sum(sapply(1:length(l), function(x) l[[x]]$sharedHits > 0))
  out.df <- data.frame("Factor" = paste0("F", 1:ncol(factorization)), "SharedTissues" = sapply(l, function(x) x$sharedHits), "pval(100)" = sapply(l, function(x) x$pval))
  write_tsv(out.df, file = paste0(args$output, "/factor_true_tissue_overlaps.ALL.txt"))
  print("Here okay...")
  if(FALSE) #me being really lazy here...
  {
    message("Calculating the AUPR per factor (weighted average)....")
    #across each factor...
    res <- NULL
    means <- c()
    for(i in 1:ncol(factorization))
    {
      r <- factorSpecificAUPR(i, trait.studies, factorization, factorization.enrichment, ldsc.reference, answer_thres = 0.05)
      n <- rep(paste0("F", i), length(r$aupr))
      ret <- cbind(n, as.numeric(r$aupr), as.numeric(r$weights))
      res <- rbind(res, ret)
      means <- c(means, r$wmean)
    }
    #compile the results into tables 
    rest <- data.frame(res)
    names(rest) <- c("factor", "aupr", "weight")
    rest$aupr <- as.numeric(as.character(rest$aupr))
    rest$weight <- as.numeric(as.character(rest$weight))
    rest$factor <- factor(rest$factor, levels = paste0("F", 1:ncol(factorization)))
    means <- data.frame("factor" = paste0("F", 1:ncol(factorization)), "wmean" = means) %>% add_row("factor" = "Average", "wmean" = mean(means,na.rm = TRUE))
    #plot
    p <- ggplot(rest, aes(x = factor, y = aupr, weight = weight)) + 
      geom_boxplot() + 
      geom_point(data=means,aes(x=factor,y=wmean), color = "coral", size = 3,shape = 4, inherit.aes=FALSE) + 
      theme_minimal(18) + xlab("Factor") + ylab("Weighted AUPR")
    
    ggsave(filename = paste0(args$output, "/pr_boxplot_factors.png"), plot = p, width = 10)
    #also a basic plot of the means
    png(paste0(args$output, "/pr_weightedmeans_factors.png") )
    barplot(names.arg = unlist(means$factor), height = means$wmean,col = "skyblue", xlab = "Factor", ylab = "weighted mean")
    dev.off() 
    #save the weighted means
    #means <- data.frame("factor" = paste0("F", 1:ncol(factorization)), "wmean" = means)
    write_tsv(x = means, file = paste0(args$output, "/pr_weightedmeans_factors.txt"))
    
    
    #this looks at overlapping enrichments between factors, and counts those tissues as true
    message(paste0("Getting the AUPR by overlap, overlap is ", args$n_overlap))
    factor.aupr <- processFactorAUPRByOverlap(factorization, trait.studies, factorization.enrichment, ldsc.reference, n_include = args$n_overlap)
    out.df <- data.frame("Factor" = paste0("F", 1:ncol(factorization)), "AUPR" = sapply(1:ncol(factorization), function(x) factor.aupr[[x]]$pr$auc.integral))
    write_tsv(x = out.df, file = paste0(args$output, "/factor_true_tissue_overlaps.", args$n_overlap, ".txt"))
  }
}

if(args$overlap_test)
{
  #this looks at overlapping enrichments between factors, and counts those tissues as true
  message(paste0("Getting the AUPR by overlap, overlap is ", args$n_overlap))
  factor.aupr <- processFactorAUPRByOverlap(factorization, trait.studies, factorization.enrichment, ldsc.reference, n_include = args$n_overlap)
  out.df <- data.frame("Factor" = paste0("F", 1:ncol(factorization)), "AUPR" = sapply(1:ncol(factorization), function(x) factor.aupr[[x]]$pr$auc.integral))
  write_tsv(x = out.df, file = paste0(args$output, "/factor_true_tissue_overlaps.", args$n_overlap, ".txt"))
}

################# Simplified metrics
if(args$all || args$simple)
{
  
  sharing <- tissueSharingMatrix(args$fdr.thresh, trait.studies,ldsc.reference)
  check.list <- apply(factorization, 2, function(x) which(abs(x) > sparsity.threshold))
  fact.scores <- sapply(check.list, function(v) simpleScore(v, sharing))
  fact.score.scaling <- sapply(check.list, function(v) choose(length(v), 2))
  fact.dist <- sapply(check.list, function(v) simpleNorm(v, sharing))
  out.df <- data.frame("Factor" = paste0("F", 1:ncol(factorization)), "Score" = unlist(fact.scores), "Dist" = unlist(fact.dist), "ScaledScore" = unlist(fact.scores)/unlist(fact.score.scaling), "ScaledAlt" = unlist(fact.scores)/sapply(check.list, length))
  out.df$Score_pval <- calcSimplePval(10000, check.list, trait.studies, sharing)
  out.df$Factor_sparsity <- apply(factorization, 2, function(x) sum(abs(x) <=  sparsity.threshold) / length(x)) #note this is confusing, counts number of about-zero entries, so higher --> more sparsity
  
  write_tsv(x = out.df, file = paste0(args$output, "/factor_simple_scores.txt"))
  #I am suspicious that denser factors are disfavored in having low p-values- like the null isn't fair.
  #plot(apply(factorization, 2, function(x) sum(x != 0)), pval)
  #check the score
  #repeat that a bunch of times
  #get a p-value
  
}
#some debugging things.....
#suspects <- c("30020_irnt", "30070_irnt", "30010_irnt")
#spy <- ldsc.reference %>%filter(Source %in% suspects) %>% filter(within_trait_FDR < 0.05)
  