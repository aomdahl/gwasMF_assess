####### Meta-analyze the factorization metrics
####### Made by Ashton Omdahl, Battle Lab, Sept 2021
rm(list = ls()) #clear the environment
pacman::p_load(data.table, tidyr, dplyr, readr, ggplot2, stringr, Xmisc, PRROC, ROCR,svMisc, stats, cowplot, ggsignif)
#######functions

if(FALSE)
{
  args <- list()
  args$factorization_handles <- "/work-zfs/abattle4/ashton/snp_networks/scratch/ldsc_all_traits/factorization_run_lists/7_k_runlist.txt"
  args$results_dir <- "/work-zfs/abattle4/ashton/snp_networks/scratch/ldsc_all_traits/results/"
  args$results_dir <- "/work-zfs/abattle4/ashton/snp_networks/gwas_decomp_ldsc/results/seed2_thresh0.9_h2-0.1_vars1e-5/factorization/K7_sept17/"
}

hackyAssignMethod <- function(handle)
{
  if(grepl("PMA", handle)){return("PMA")
  } else if(grepl("sPCA", handle) | grepl("sparsePCA", handle)) {return("sparsePCA")
  } else if(grepl("ssvd", handle) || grepl("sSVD", handle)) {return("ssvd")
  } else if(grepl("flash", handle)) {return("flashR")
  } else if(grepl("B_SE", handle) | grepl("beta_se", handle)) {return("gwasMF")
  } else if(grepl("backfit", handle)) {return("flashR")
  } else if(grepl("ashr", handle)) {return("flashR")
  } else if(grepl("PCA", handle)) {return("PCA")
  } else {return("unknown")}
}
################### MAIN ###################

parser <- ArgumentParser$new()
parser$add_description("Script to examine results")
parser$add_argument("--factorization_handles", type = 'character', help = "File containing factorization handles" )
parser$add_argument("--results_dir", type = 'character', help = "path to results dir" )
parser$add_argument("--test_type", type = 'character', help = "Specify the kind of assessment you are doing- general, or celltype", default = "celltype" )
parser$add_argument("--output", type = 'character', help = "output directory to write to")
parser$add_argument('--help',type='logical',action='store_true',help='Print the help page')

parser$helpme()
args <- parser$get_args()
handles <- scan(args$factorization_handles, what = character())
paths <- paste0(args$results_dir, "/", handles, "/")


if(args$test_type == "celltype")
{
  handles <- scan(args$factorization_handles, what = character())
  #Overlap measures
  n_overlaps <- do.call("rbind", lapply(handles, function(f) fread(paste0(paste0(args$results_dir, f), "/factor_true_tissue_overlaps.2.txt")) %>% mutate("Factorization" = f))) %>% mutate("method" = factor(sapply(.$Factorization, hackyAssignMethod))) %>% arrange(method)
  
  #with.nas <- data.frame("Factorization" = unique((n_overlaps %>% filter(is.na(AUPR)))$Factorization), "AUPR" = 1)
  
  ggplot(data = n_overlaps, aes(x = reorder(Factorization, as.integer(method)), y = AUPR, fill = method)) + geom_boxplot() + theme_minimal(17) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + xlab("AUPR across factors") #+ 
  geom_text(data = with.nas, aes(label = "*", y = AUPR, x = Factorization),inherit.aes=FALSE, color = "red")
  
  #simple measures
  simple <- do.call("rbind", lapply(handles, function(f) fread(paste0(paste0(args$results_dir, f), "/factor_simple_scores.txt")) %>% mutate("Factorization" = f))) %>% mutate("method" = factor(sapply(.$Factorization, hackyAssignMethod))) %>% arrange(method)
  simple$Factor <- factor(simple$Factor, levels = paste0("F", 1:20))
  
  simple.means <- simple %>% group_by(Factorization, method) %>% summarize("avg" = mean(Score))
  ggplot(data = simple.means, aes(x = Factorization, y = avg, fill = method)) + geom_bar(stat = "identity") +
    theme_minimal(15) + ggtitle("Avg score") + theme(axis.text.x = element_text(angle = 90))

  #overall performance
  sums.overall <- simple %>% group_by(Factorization, method) %>% summarize("tot" = sum(Score))
  ggplot(data = sums.overall, aes(x = Factorization, y = tot, fill = method)) + geom_bar(stat = "identity") + 
    theme_minimal(15) + ggtitle("Overall score (summed)") + theme(axis.text.x = element_text(angle = 90))
  
  #make a boxplot YOU NINNY
  ggplot(data = simple, aes(x = reorder(Factorization, as.integer(method)), y = Score, fill = method)) + geom_boxplot() + theme_minimal(15) + ggtitle("Scores across factors by method") + theme(axis.text.x = element_text(angle = 90)) + xlab("Factorization method")
  
  #boxplot of pvalues
  ggplot(data = simple, aes(x = reorder(Factorization, as.integer(method)), y = Score_pval, fill = method)) + geom_boxplot() + theme_minimal(15) + ggtitle("Score Pvals across factors by method") + theme(axis.text.x = element_text(angle = 90)) + xlab("Factorization method")
  
  #plot of pvals vs density
  res <- lm(Score_pval ~ Factor_sparsity, data = simple)
  ggplot(data = simple, aes(x = Factor_sparsity, y = Score_pval)) + geom_point() +geom_smooth(method = "lm") +
    theme_minimal(15) + annotate(geom = "text",x = 0.75, y = 0.5, label = paste0("p = ", round(summary(res)$coefficients[8],digits = 3))) + geom_hline(yintercept = 0.05, color = "red", linetype = 'dashed')
  hist(simple$Score_pval, main = "Distribution of p-values", xlab = "Score pvalue", breaks = 30)
}
if(args$test_type == "general")
{
  #Frob norm
  fs <- list.files(args$results_dir, pattern = "*recon_error.txt")
  n_frob <- do.call("rbind", lapply(fs, function(f) fread(paste0(args$results_dir,"/", f)) %>% mutate("Method" = str_split(f, pattern = "\\.")[[1]][1]) %>% mutate("Type" = hackyAssignMethod(Method))))
  str_split(fs, pattern = "\\.")
  ggplot(n_frob, aes(x = Method, y = Frobenius_norm, fill = Type)) + geom_bar(stat = "identity") + theme_minimal(15) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ggtitle("Reconstruction error (frobenius norm)")
  ggsave(paste0(args$output, "/frob.png"))
  #sparsity and pve
  sp.pve <- list.files(args$results_dir, pattern = "*pve.txt")
  n_sp <- do.call("rbind", c(lapply(sp.pve, function(f) fread(paste0(args$results_dir,"/", f)) %>% mutate("Method" = str_split(f, pattern = "\\.")[[1]][1])), fill = TRUE))
  str_split(fs, pattern = "\\.")
  ggplot(n_sp, aes(x = Factor, y = PVE, color = Method)) + geom_line(aes(group = Method)) + theme_minimal(15) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ggtitle("PVE (estimates)")
  ggsave(paste0(args$output, "/pve.png"))
  
  ggplot(n_sp, aes(x = Factor, y = Sparsity, color = Method)) + geom_line(aes(group = Method)) + theme_minimal(15) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ggtitle("Sparsity (|x| < 1e=4)")
}
ggsave(paste0(args$output, "/sparsity.png"))

