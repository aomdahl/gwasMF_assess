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