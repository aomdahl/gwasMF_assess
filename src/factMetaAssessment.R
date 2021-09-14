####### Meta-analyze the factorization metrics
####### Made by Ashton Omdahl, Battle Lab, Sept 2021
rm(list = ls()) #clear the environment
install.packages("ggsignif")
pacman::p_load(data.table, tidyr, dplyr, readr, ggplot2, stringr, Xmisc, PRROC, ROCR,svMisc, stats, cowplot, ggsignif)
#######functions

if(FALSE)
{
  args <- list()
  args$factorization_handles <- "/work-zfs/abattle4/ashton/snp_networks/scratch/ldsc_all_traits/factorization_run_lists/sept2.txt"
  args$results_dir <- "/work-zfs/abattle4/ashton/snp_networks/scratch/ldsc_all_traits/results/"
}

hackyAssignMethod <- function(handle)
{
  if(grepl("PMA", handle)){return("PMA")
  } else if(grepl("sPCA", handle) | grepl("sparsePCA", handle)) {return("sparsePCA")
  } else if(grepl("ssvd", handle)) {return("ssvd")
  } else if(grepl("flash", handle)) {return("flashR")
  } else if(grepl("B_SE", handle) | grepl("beta_se", handle)) {return("gwasMF")
  } else if(grepl("backfit", handle)) {return("flashR")
  } else if(grepl("ashr", handle)) {return("flashR")
  } else {return("unknown")}
}
################### MAIN ###################

parser <- ArgumentParser$new()
parser$add_description("Script to examine results")
parser$add_argument("--factorization_handles", type = 'character', help = "File containing factorization handles" )
parser$add_argument("--results_dir", type = 'character', help = "path to results dir" )
parser$add_argument("--output", type = 'character', help = "output directory to write to")

parser$helpme()
args <- parser$get_args()

handles <- scan(args$factorization_handles, what = character())
#Overlap measures
n_overlaps <- do.call("rbind", lapply(handles, function(f) fread(paste0(paste0(args$results_dir, f), "/factor_true_tissue_overlaps.2.txt")) %>% mutate("Factorization" = f))) %>% mutate("method" = factor(sapply(.$Factorization, hackyAssignMethod))) %>% arrange(method)

with.nas <- data.frame("Factorization" = unique((n_overlaps %>% filter(is.na(AUPR)))$Factorization), "AUPR" = 1)

ggplot(data = n_overlaps, aes(x = reorder(Factorization, as.integer(method)), y = AUPR, fill = method)) + geom_boxplot() + theme_minimal(17) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + xlab("AUPR across factors") + 
  geom_text(data = with.nas, aes(label = "*", y = AUPR, x = Factorization),inherit.aes=FALSE, color = "red")

#simple measures
simple <- do.call("rbind", lapply(handles, function(f) fread(paste0(paste0(args$results_dir, f), "/factor_simple_scores.txt")) %>% mutate("Factorization" = f))) %>% mutate("method" = factor(sapply(.$Factorization, hackyAssignMethod))) %>% arrange(method)
simple$Factor <- factor(simple$Factor, levels = paste0("F", 1:20))

#p.flash <- ggplot(data = simple %>% filter(method == "flashR"), aes(x = Factor, y = Score, color = method)) + geom_line(aes#(group = Factorization)) + theme_minimal(15)
#p.svd <- ggplot(data = simple%>% filter(method %in% c("ssvd", "sparsePCA")), aes(x = Factor, y = Score, color = method)) + geom_line(aes(group = Factorization)) + theme_minimal(15)
#p.pma <- ggplot(data = simple %>% filter(method == "PMA"), aes(x = Factor, y = Score, color = method)) + geom_line(aes(group = Factorization)) + theme_minimal(15)
#p.gwasmf <- ggplot(data = simple %>% filter(method == "gwasMF"), aes(x = Factor, y = Score, color = method)) + geom_line(aes(group = Factorization)) + theme_minimal(15)
#plot_grid(p.flash, p.svd, p.pma, p.gwasmf)
#ggplot(data = simple, aes(x = Factor, y = ScaledAlt, color = method)) + geom_line(aes(group = Factorization)) + theme_minimal(15) + ggtitle("Scaled by # traits")

simple.means <- simple %>% group_by(Factorization, method) %>% summarize("avg" = mean(Score))
ggplot(data = simple.means, aes(x = Factorization, y = avg, fill = method)) + geom_bar(stat = "identity") +
  theme_minimal(15) + ggtitle("Avg score") + theme(axis.text.x = element_text(angle = 90))

#means.2 <- simple %>% group_by(Factorization, method) %>% summarize("avg" = mean(ScaledAlt))
#ggplot(data = means.2, aes(x = Factorization, y = avg, fill = method)) + geom_bar(stat = "identity") + theme_minimal(15) + ggtitle("Mean scaled by # traits")
#overall performance
sums.overall <- simple %>% group_by(Factorization, method) %>% summarize("tot" = sum(Score))
ggplot(data = sums.overall, aes(x = Factorization, y = tot, fill = method)) + geom_bar(stat = "identity") + 
  theme_minimal(15) + ggtitle("Overall score (summed)") + theme(axis.text.x = element_text(angle = 90))

#TODO- sum over all facctors somehow?
