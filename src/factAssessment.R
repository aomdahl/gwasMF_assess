####### Assess the quality of a factorization
####### Made by Ashton Omdahl, Battle Lab, Aug 2021
rm(list = ls()) #clear the environment
pacman::p_load(data.table, tidyr, dplyr, readr, ggplot2, stringr, Xmisc, PRROC, ROCR,svMisc, stats, cowplot)
#######functions

#plot the pr curve using existing package
#@param trait.id: identifier of the trait you want to lookup
#@param factorization: the actual factorization. Used to select which factors are relevant
#@param factorization.enrichment: the data listing 
areaPRPlot <- function(ind, trait.studies, factorization, factorization.enrichment, ldsc.reference, answer_thres = 0.05)
{
  trait.id <- trait.studies[ind]
  #Gets all the info for the traits
  fact.enrichment <- filterFactorization(trait.index = ind, fact.matrix = factorization, fact.enrichment = factorization.enrichment, FDR = 1, overall_FDR = TRUE) #give me all the hits
  
  #get the correct answers
  all.tissues <- unique((ldsc.reference %>% arrange(tissue))$tissue)
  ref.enrichment <- filterReference(ldsc.reference, trait.id, FDR = answer_thres, overall_FDR = FALSE) 
  
  #Check for some weird edge cases...
  if(nrow(fact.enrichment) == 0) #we have an empty factor, no enrichments in the factor
  {
    #est <- data.frame("names" = all.tissues, "labels" = ref.panel, "probs" = runif(length(all.tissues), 0, 0.05) )
    psuedo <- list()
    psuedo$auc.integral <- 0
    psuedo$auc.davis.goadrich <- 0
    return(psuedo)
    
  }else{ #we do have something..
    answer.fact <-  fact.enrichment %>% group_by(tissue) %>% slice(which.min(overall_FDR)) %>% mutate("class_prob" = 1-overall_FDR) %>% ungroup() %>% arrange(tissue)
  }
  
  #another bunch of edge cases, this time in ref enrichment
  if(!any(all.tissues %in% ref.enrichment$tissue)) #this trait had no enrichments in LDSC
  {
    ref.panel <- rep(0, length(all.tissues)) #this is the case when its NA....
    est <- data.frame("names" = all.tissues, "labels" = ref.panel, "probs" = answer.fact$class_prob ) %>% add_row(names = "extra", labels = 1, probs = 0.5)
  }else {
    ref.panel <- as.integer(all.tissues %in% ref.enrichment$tissue)
  est <- data.frame("names" = all.tissues, "labels" = ref.panel, "probs" = answer.fact$class_prob )  #TODO: evaluate this. is this correct?
  }

  #full list of categoryie
  back <- PRROC::pr.curve(scores.class0=est[est$labels==1,]$probs,
                          scores.class1=est[est$labels==0,]$probs,
                          curve=T)
  return(back)
}


#Read in the LDSC output directory
#currently being lazy and just assuming reference same for all...
readInLDSC <- function(dir, ext = ".cell_type_results.txt", type = "ldsc reference")
{
  file.list <- list.files(dir, pattern = ext)
  
  if(type == "ldsc reference")
  {
    all <- lapply(file.list, function(x) fread(paste0(dir, x)) %>% arrange("Name") %>% mutate("Source" = gsub(x = x,pattern = ".cell_type_results.txt", replacement = "")))
  }else
  {
    all <- lapply(file.list, function(x) fread(paste0(dir, x)) %>% arrange("Name") %>% mutate("Source" = str_extract(str_split(x, pattern = "\\.")[[1]][1], pattern = "F\\d+"))) 
  }
  todos <- plyr::rbind.fill(all) %>% mutate("-log10(P)" =-log10(Coefficient_P_value))
  labels <- fread("/work-zfs/abattle4/ashton/snp_networks/custom_l1_factorization/ldsc_reference/finucane_2018_supp_table7.ashton.csv")
  #Kyped the following code from Rebecca Keener on 6/2/2021
  yes<-subset(labels, labels$entex=="Yes")
  yes$Name<-paste(yes$tissue, "_ENTEX__", yes$mark, sep="")
  no<-subset(labels, labels$entex=="No")
  no$Name<-paste(no$tissue, no$mark, sep="__")
  labels<-rbind(yes, no)
  substitution_regex = "__"
  substitution_regex_to = "-"
  
  ret <- left_join(todos, labels, by  = "Name") %>% arrange(new_category) %>% select(-entex) %>% group_by(Source) %>% 
    mutate("within_trait_FDR" = p.adjust(Coefficient_P_value, method = "fdr")) %>% ungroup()
  ret$overall_FDR <- p.adjust(ret$Coefficient_P_value, method = "fdr")
  ret
}
#Filter the reference data from LDSC
#@param ldsc.ref: the read in ldsc reference information
#@param trait.id: identifier for the trait
#@param FDR: the FDR of what we accept as an answer
#@param overall_FDR: if we use the trait-specific FDR or overall FDR. Trait specific is default, since examining only one trait at a time.
filterReference <- function(ldsc.ref, trait.id, FDR = 0.01, overall_FDR = FALSE)
{
  n <-  ldsc.ref %>% filter(Source == trait.id)
  if(overall_FDR)
  {
    n %>% ungroup() %>% filter(overall_FDR < FDR) %>% select(tissue, mark, new_category)
    
  } else{
    n %>% ungroup() %>% filter(within_trait_FDR < FDR) %>% select(tissue, mark, new_category)
  }
}

#Filter the LDSC analysis of the factorization
#@param trait.index: the index of the trait of interest
#@param fact.matrix: the matrix of factorization information
#@param FDR: fact.enrichmentL: the table of LDSC enrichment for it.
#@param overall_FDR: if we use the trait-specific FDR or overall FDR. Overall is default, since we look at all factors at once (?)
filterFactorization <- function(trait.index, fact.matrix, fact.enrichment, FDR = 0.05, overall_FDR = TRUE)
{
  factors <- which(fact.matrix[trait.index,] != 0)
  if(length(factors) == 0)
  {
    #no factors significant here; return an empty table
    print(paste0("No factors significant for ", trait.index))
    return(fact.enrichment %>% filter(overall_FDR < 0) %>% select(Source, tissue, mark, new_category))
  }
  relative_weights <- fact.matrix[trait.index,..factors]^2/sum(fact.matrix[trait.index,..factors]^2)
  filtered.fact <- fact.enrichment %>% filter(Source %in% paste0("F", factors))
  if(overall_FDR)
  {
    filtered.fact <- filtered.fact %>% filter(overall_FDR < FDR) %>% select(Source, tissue, mark, new_category, overall_FDR)
  }else
  {
    filtered.fact <- filtered.fact %>% filter(within_trait_fdr < FDR) %>% select(Source, tissue, mark, new_category, overall_FDR)
  }
  return(filtered.fact)
}

checkTraitEnrichments <- function(trait.list,trait.ids,  ldsc.reference, answer_thresh)
{
  relevant.traits <- c()
  #count which ones have real enrichments at the answer threshold
  for(t in trait.list)
  {
    trait.id <- trait.studies[t]
    #Identify all the "true" enrichment for said trait
    ref.enrichment <- filterReference(ldsc.reference, trait.id, FDR = answer_thresh, overall_FDR = FALSE) %>% arrange(tissue)
    if(nrow(ref.enrichment) != 0)
    {
      relevant.traits <- c(relevant.traits, t)
    }else
    {
      print("omitting trait #....")
      print(t)
    }
  }
  return(relevant.traits)
}

#@PARAM f.ind: the index of the factor we are looking at
#@param trait.studies: the list of trait study identifiers
#@param factorization: the actual factorization matrix
#@param factorization.enrichment: the enrichment information for each factor
#@param ldsc.reference: the ldsc enrichment information for each trait (our answers/reference
#@param answer_thresh: the threshold at which enrichments are counted as true
#Note that here we depart from previous method by accounting for within-"trait" FDR (i.e. across a single factor), not across all factors.
#The reason for this is that we are interested only in the performance of a single factor.
#This can easily be changed, see "factor.enrichments <- ..."
factorSpecificAUPR <- function(f.ind, trait.studies, factorization, factorization.enrichment, ldsc.reference, answer_thres = 0.05)
{
  pr.list <- list()
  factor.id <- paste0("F", f.ind)
  #get the traits in the factor that are nonzero
  relevant.traits <- which(factorization[,..f.ind] != 0)
  #determine all of the tissues which are enriched for that factor and their FDR score. Pick only the highest marker per tissue.
  factor.enrichments <- factorization.enrichment %>% filter(Source == factor.id) %>% mutate("class_prob" = 1-within_trait_FDR) %>% group_by(tissue) %>% slice(which.max(class_prob)) %>% 
    ungroup() %>% select(tissue, mark, new_category, class_prob) %>% arrange(tissue) #doing here for within_trait FDR? *** very important
  all.tissues <- unique((factor.enrichments %>% arrange(tissue))$tissue)
  #FOR EACH of those traits
  
  #print(relevant.traits)
  #relevant.traits <- checkTraitEnrichments(relevant.traits,trait.ids,  ldsc.reference, answer_thres)
  #print(relevant.traits)
  #readline()
  
  for(t in relevant.traits)
  {
    trait.id <- trait.studies[t]
    #Identify all the "true" enrichment for said trait
    ref.enrichment <- filterReference(ldsc.reference, trait.id, FDR = answer_thres, overall_FDR = FALSE) %>% arrange(tissue)
    #Of those, which are enriched in Factor matrix?
    ref.panel <- as.integer(all.tissues %in% unique(ref.enrichment$tissue))
    #the label is 1 if a trait is enriched for that tissue at FDR = answer_thres, 0 otherwise
    est <- data.frame("names" = all.tissues, "labels" = ref.panel, "probs" = factor.enrichments$class_prob )
    
    #get the area under the PR curve
    #Aug 27 executive call:
      #if a trait has no enrichments at all, then we don't include it.
    pr.list[[t]] <- PRROC::pr.curve(scores.class0=est[est$labels==1,]$probs,
                                    scores.class1=est[est$labels==0,]$probs,
                                    curve=T)
    #you know, maybe the play is for tissues with no enrichments at all, throw those out.
    #no, we want to be able to detect
    if(FALSE)
    #if(is.na(pr.list[[t]]$auc.integral)) #this means there were no enrichments in the true tissue....
    {
      print("Found na")
      print(paste0("factor", f.ind))
      print(paste0("Trait", t))
      print(est)
      readline()
    }
  }
  #calculate a weighted average...
  library(stats)
  pr.scores <- unlist(lapply(relevant.traits, function(i) pr.list[[i]]$auc.integral)) #extract the actual PR scores we've calculated
  drop <- which(is.na(pr.scores)) #we don't include them in our mean calculation.
  if(length(drop) < 1)
  {
    weights <-  unlist(factorization[relevant.traits,..f.ind]^2 / sum(factorization[relevant.traits,..f.ind]^2)) #extract the weights
  }else
  {
    weights <-  unlist(factorization[relevant.traits,..f.ind]^2 / sum(factorization[relevant.traits,..f.ind]^2))[-drop] #extract the weights
    pr.scores <- pr.scores[-drop]
  }

  back <- list()
  back$traits <- relevant.traits
  back$aupr <- pr.scores
  back$weights <- weights
  #t <- weighted.mean(x = pr.scores, w = weights)
  back$wmean <- weighted.mean(x = pr.scores, w = weights)
  
  return(back)
}

#@PARAM f.ind: the index of the factor we are looking at
#@param trait.studies: the list of trait study identifiers
#@param factorization: the actual factorization matrix
#@param factorization.enrichment: the enrichment information for each factor
#@param ldsc.reference: the ldsc enrichment information for each trait (our answers/reference
#@param answer_thresh: the threshold at which enrichments are counted as true
#Note that here we depart from previous method by accounting for within-"trait" FDR (i.e. across a single factor), not across all factors.
#This differs from the above method in that we look at perf
factorSpecificAUPRShared <- function(f.ind, trait.studies, factorization, factorization.enrichment, ldsc.reference, answer_thres = 0.1)
{
  pr.list <- list()
  factor.id <- paste0("F", f.ind)
  #get the traits in the factor that are nonzero
  relevant.traits <- which(factorization[,..f.ind] != 0)
  #determine all of the tissues which are enriched for that factor and their FDR score. Pick only the highest marker per tissue.
  factor.enrichments <- factorization.enrichment %>% filter(Source == factor.id) %>% mutate("class_prob" = 1-within_trait_FDR) %>% group_by(tissue) %>% slice(which.max(class_prob)) %>% 
    ungroup() %>% select(tissue, mark, new_category, class_prob) %>% arrange(tissue) #doing here for within_trait FDR? *** very important
  
  all.tissues <- unique((factor.enrichments %>% arrange(tissue))$tissue)
  
  #FOR EACH of those traits, track which enrichments shared by all (intersect)
  for(i in 1:length(relevant.traits)){
    t = relevant.traits[i]
    trait.id <- trait.studies[t]
    ref.enrichment <- filterReference(ldsc.reference, trait.id, FDR = answer_thres, overall_FDR = FALSE) %>% arrange(tissue)
    if(i == 1)
    {
      tissue.relevant.traits <- unique(ref.enrichment$tissue)
      #print(tissue.relevant.traits)
    } else  {
      tissue.relevant.traits <- intersect(tissue.relevant.traits,unique(ref.enrichment$tissue))
      #print(tissue.relevant.traits)
    }
  }
  ref.panel <- as.integer(all.tissues %in% tissue.relevant.traits)
  #the label is 1 if a trait is enriched for that tissue at FDR = answer_thres, 0 otherwise
  est <- data.frame("names" = all.tissues, "labels" = ref.panel, 
                    "probs" = factor.enrichments$class_prob )
  
  #get the area under the PR curve
  pr.list <- PRROC::pr.curve(scores.class0=est[est$labels==1,]$probs,
                             scores.class1=est[est$labels==0,]$probs,
                             curve=T)
  back <- list()
  back$FactorTraits <- relevant.traits
  back$aupr <- pr.list
  back$sharedHits <- sum(est$labels)
  back$sharedTissues <- tissue.relevant.traits
  return(back)
}

#HELPER function: For each tissue, update the number of traits that contain it.
#Helper for @function factorTraitSharingMatrix
#@param sharing.dat: sharing data
updateSharingDat <- function(sharing.dat, tissue, curr.list, trait,id)
{
  for(t in tissue)
  {
    if(t %in% curr.list)
    {
      sharing.dat[[t]]$count = sharing.dat[[t]]$count + 1
      sharing.dat[[t]]$traits <- paste0(sharing.dat[[t]]$traits,",", trait)
      sharing.dat[[t]]$trait.ids <- paste0(sharing.dat[[t]]$trait.ids,",", id)
    }
  }
  return(sharing.dat)
}

#Count the tissue sharing at a specific threshold for a specific factor
#@return ret$matrix: matrix containing the pairwise count of overlapping tissues between traits in a factor
#@return ret$tissue.df: data frame containing the number of overlaps for each tissue in the factor, and which traits it appears in.
#@param f.ind: factor index
#@param trait.studies: trait names
#@param factorization: the factorization F matrix
#@param factorization.enrichment: the actual LDSc enrichment of the factorization, one table
#@param ldsc.refernece: the true known enrichemtns for each trait
#@param answer_thresh: FDR threshold (PER TRAIT) for true ldsc reference.
factorTraitSharingMatrix <- function(f.ind, trait.studies, factorization, factorization.enrichment, ldsc.reference, answer_thres = 0.05)
{
  factor.id <- paste0("F", f.ind)
  relevant.traits <- which(factorization[,..f.ind] != 0)
  all.tissues <- unique((factorization.enrichment %>% arrange(tissue))$tissue)
  
  #for n > 1:
  #join the list of all of themat the threshold....
  tissue.relevant.traits <- NULL
  shared <- matrix(NA, length(relevant.traits),length(relevant.traits))
  sharingDat <- list()
  #initialize list
  for(t in all.tissues)
  {
    sharingDat[[t]] <- list()
    sharingDat[[t]]$count <- 0
    sharingDat[[t]]$traits <- ""
    sharingDat[[t]]$trait.ids <- ""
  }
  #For all the relevant traits for a given factor
  for(i in 1:length(relevant.traits)){
    t = relevant.traits[i]
    trait.id <- trait.studies[t]
    trait.name <- trait.names[t]
    ref.enrichment <- filterReference(ldsc.reference, trait.id, FDR = answer_thres, overall_FDR = FALSE) %>% arrange(tissue)
    for (j in relevant.traits[-i])
    {
      o <- trait.studies[j]
      alt.index = which(relevant.traits == j)
      alt.ref.enrichment <- filterReference(ldsc.reference, o, FDR = answer_thres, overall_FDR = FALSE) %>% arrange(tissue)
      #shared overlap
      shared[i,alt.index] <- length(intersect(unique(alt.ref.enrichment$tissue), unique(ref.enrichment$tissue)))
    }
    sharingDat <- updateSharingDat(sharingDat, all.tissues, unique(ref.enrichment$tissue), trait.name, trait.id)
  }
  #Finally, clean up the sharing dat into a nice table
  df <- data.frame(do.call("rbind", sharingDat))
  df$count <- as.integer(df$count)
  df$tissue <- rownames(df)
  df$traits <- gsub(x = df$traits,pattern = "^,", replacement = "")
  df$trait.ids <- gsub(x =   df$trait.ids ,pattern = "^,", replacement = "")
  ret_list <- list(); ret_list$matrix <- shared; ret_list$tissue.df <- df %>% arrange(-count)
  return(ret_list)
  
}

#Function to calculate the AUPR for each factor based on overlapping tissue enrichments between studies
#@preturn ret: list for indices 1:K containing the pr curve information for each factor, and with l$df containing the tissue enrichemnts in tabular format.
#@param n_include: the requirement for overlaps to count a tissue. Default is 2 (that is, if a tissue is enriched pariwise in 2 traits, we count that as "true")
#@param trait.studies: trait names
#@param factorization: the factorization F matrix
#@param factorization.enrichment: the actual LDSc enrichment of the factorization, one table
#@param ldsc.refernece: the true known enrichemtns for each trait
#@param answer_thresh: FDR threshold (PER TRAIT) for true ldsc reference.
processFactorAUPRByOverlap <- function(factorization, trait.studies, factorization.enrichment, ldsc.reference, n_include = 2){
  ret <- list()
  all_f <- NULL
  for(f in 1:ncol(factorization))
  {
    all.tissues <-unique((ldsc.reference %>% arrange(tissue))$tissue)

    ret[[f]] <- list()
    r <-  factorTraitSharingMatrix(f, trait.studies, factorization, factorization.enrichment, ldsc.reference)
    
    tab <- r$tissue.df %>% filter(count >= n_include) %>% mutate("Factor" = paste0("F", f), "numTraits" = length(which(factorization[,..f] != 0)))
    ref.panel <- as.integer(all.tissues %in% tab$tissue)
    factor.enrichments <- factorization.enrichment %>% filter(Source == paste0("F", f)) %>% mutate("class_prob" = 1-within_trait_FDR) %>% group_by(tissue) %>% slice(which.max(class_prob)) %>% 
      ungroup() %>% select(tissue, mark, new_category, class_prob) %>% arrange(tissue) #doing here for within_trait FDR? *** very important
    est <- data.frame("names" = all.tissues, "labels" = ref.panel, 
                      "probs" = factor.enrichments$class_prob )

    #get the area under the PR curve
    pr.list <- PRROC::pr.curve(scores.class0=est[est$labels==1,]$probs,
                               scores.class1=est[est$labels==0,]$probs,
                               curve=TRUE)
    if(nrow(r$tissue.df %>% filter(count >= n_include)) < 1)
    {
      print(paste0("For current factorization, factor ", f," has no overlapping enrichments. Setting AUPR to 0"))
      #In this case, make the PRcurve 0 
      #print(r$tissue.df)
      #readline()
      pr.list$auc.integral <- 0
    }
    if(is.na(pr.list$auc.integral))
    {
      print(paste0("For current factorization, factor ", f," we have an AUPR of NA"))
      print(est)
      readline()
    }
    
    all_f <- rbind(all_f, r$tissue.df %>% filter(count > 1) %>% mutate("Factor" = paste0("F", f), "numTraits" = length(which(factorization[,..f] != 0))))
    ret[[f]]$pr <- pr.list
  }
  ret$df <- all_f
  return(ret)
}
 

###Simple helpers
#Create a binary matrix indicating if traits share some tissue at a set fdr.thresh
tissueSharingMatrix <- function(fdr.thresh, trait.studies, ldsc.reference)
{
  sharing <- matrix(0, length(trait.studies), length(trait.studies))
  for(t in trait.studies)
  {
    curr.t <- unique((ldsc.reference %>% filter(Source == t, within_trait_FDR < fdr.thresh))$tissue)
    if(length(curr.t) == 0) #no matches can occur
    {
      next
    }
    matches <-ldsc.reference %>% filter(within_trait_FDR < fdr.thresh, tissue %in% curr.t) %>% group_by(Source) %>% 
      slice(which.min(within_trait_FDR)) %>% ungroup()
    #^ we allow for matchin with onese self, makes downstream easier. Just omit it.
    curr.i <- which(t == trait.studies)
    match.i <- sapply(matches$Source, function(j) which(j == trait.studies))
    sharing[curr.i, match.i] <- 1
  }
  return(sharing)
}
#@param vect: the indices of the traits in the group
#@param m: the matrix of pairwise tissue sharing across traits
simpleScore <- function(vect, m)
{
  true <- m[vect,vect] #what the actual relationship between these tissues is
  true[lower.tri(true, diag = TRUE)] <- NA #only want the upper triangle
  #assume our guess is all in here should be one
  sum(true == 1,na.rm = TRUE) - sum(true == 0, na.rm = TRUE)
  #The ones we guessed right in the group - the ones we didn't
}
#I think maybe this is the way to go?
simpleNorm <- function(vect, m)
{
  true <- m[vect,vect] #what the actual relationship between these tissues is
  t <- true[lower.tri(true, diag = FALSE)] #only want the upper triangle
  all.ones <- matrix(1, length(vect), length(vect))
  p <- all.ones[lower.tri(all.ones, diag = FALSE)]
  sqrt(sum((t-p)^2))
  #The ones we guessed right in the group vs the ones we didn't
}
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
parser$add_argument("--all", type = 'logical', default = FALSE, help = "Run all analysis metrics.")
parser$add_argument("--trait_specific", type = 'logical', default = FALSE, help = "Run trait_specific analysis metrics.")
parser$add_argument("--factor_specific", type = 'logical', default = FALSE, help = "Run factor-specific analysis metrics.")
parser$add_argument("--simple", type = 'logical', default = FALSE, help = "Run simple tests.")
parser$add_argument("--fdr.thresh", type = 'double', default = 0.05, help = "True value cutoff.")

parser$add_argument("--overlap_test", type = 'logical', default = FALSE, help = "Run overlap test.")
parser$add_argument('--help',type='logical',action='store_true',help='Print the help page')
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

trait.names <- scan("/work-zfs/abattle4/ashton/snp_networks/gwas_decomp_ldsc/trait_selections/seed2_thresh0.9_h2-0.1.names.tsv", what = character())
trait.studies <- scan("/work-zfs/abattle4/ashton/snp_networks/gwas_decomp_ldsc/trait_selections/seed2_thresh0.9_h2-0.1.studies.tsv", what = character())
K <- ncol(factorization)
ntraits <- nrow(factorization)
res <- matrix(NA, ntraits, 2)

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
  
  message("Calculating the AUPR per factor (weighted average)....")
  if(FALSE) #this was the old way of doing it, now have a better way..
  {
    factor.avgs <- sapply(1:K, function(i) weighted.mean(x = res[,1], w = unlist(factorization[,..i]^2)))
    png(paste0(args$output, "/pr_per_factor.png") )
    barplot(factor.avgs, names.arg = 1:K, ylab = "Factors", xlab = "Average AUPR", horiz = TRUE)
    dev.off()
    by.factor <- data.frame("Factor" = paste0("F", 1:K), "Weighted_avg" = factor.avgs)
    write_tsv(x = by.factor, file= paste0(args$output, "/per_factor_aupr.tsv"))
  }
  
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
  
  #AUPR per factor by sharing
  #Number of shared enrichments within a factor..... (how good is the grouping?)
  l <- list()
  for(i in 1:ncol(factorization))
  {
    l[[i]] <- factorSpecificAUPRShared(i, trait.studies, factorization, factorization.enrichment, ldsc.reference, answer_thres = 0.1)
  } 
  num.factors.with.enrichemnts <- sum(sapply(1:length(l), function(x) l[[x]]$sharedHits > 0))
  fileConn<-file(paste0(args$output, "/factor_true_tissue_overlaps.ALL.txt"))
  writeLines(paste0("Number of factors for which ALL grouped traits shared a detected overlapping tissue association: ", num.factors.with.enrichemnts), fileConn)
  close(fileConn)
  
    #this looks at overlapping enrichments between factors, and counts those tissues as true
  message(paste0("Getting the AUPR by overlap, overlap is ", args$n_overlap))
  factor.aupr <- processFactorAUPRByOverlap(factorization, trait.studies, factorization.enrichment, ldsc.reference, n_include = args$n_overlap)
  out.df <- data.frame("Factor" = paste0("F", 1:ncol(factorization)), "AUPR" = sapply(1:ncol(factorization), function(x) factor.aupr[[x]]$pr$auc.integral))
  write_tsv(x = out.df, file = paste0(args$output, "/factor_true_tissue_overlaps.", args$n_overlap, ".txt"))
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
  check.list <- apply(factorization, 2, function(x) which(abs(x) > 5e-3))
  fact.scores <- sapply(check.list, function(v) simpleScore(v, sharing))
  fact.score.scaling <- sapply(check.list, function(v) choose(length(v), 2))
  fact.dist <- sapply(check.list, function(v) simpleNorm(v, sharing))
  out.df <- data.frame("Factor" = paste0("F", 1:ncol(factorization)), "Score" = unlist(fact.scores), "Dist" = unlist(fact.dist), "ScaledScore" = unlist(fact.scores)/unlist(fact.score.scaling), "ScaledAlt" = unlist(fact.scores)/sapply(check.list, length))
  write_tsv(x = out.df, file = paste0(args$output, "/factor_simple_scores.txt"))
}
#some debugging things.....
#suspects <- c("30020_irnt", "30070_irnt", "30010_irnt")
#spy <- ldsc.reference %>%filter(Source %in% suspects) %>% filter(within_trait_FDR < 0.05)
  