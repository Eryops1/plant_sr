# arguments for manual parallelization on server
#args = commandArgs(trailingOnly=TRUE)
#x = as.numeric(args[1])

#setwd("~/Documents/WOLF/PROJECTS/58 World Checklist paper/analyses 2019")
#setwd("/data_vol/wolf/WCSP_paper/")

library(ape)
library(parallel)
source("plant_sr/functions.R")

# retrieve link tables, phylogenies, and wcsp data (from add_species.R)
load("plant_sr_data/MATCHES_a_a.RData")
# includes MATCHES_a_a, MATCHES_a_a_a, phylo_a, phylo_a_a_a, wcsp
## MATCHES_a_a: dataframe with columns tip | conservative
## MATCHES_a_a_a: shorter df with columns tip | conservative

load("plant_sr_data/MATCHES_b_a.RData")
# MATCHES_b_a, MATCHES_b_a_a, 
rm(MATCHES_a_a, MATCHES_b_a, phylo_a, phylo_b)

# retrieve community matrix (from process_geography.R)
load("plant_sr_data/comm.RData")

####################################################
# 1. Computing tree-level variables (RD and EDGES) #
####################################################

 # keep <- sample(1:length(phylo_a$tip.label), 300)
 # drop <- 1:length(phylo_a$tip.label)
 # drop <- drop[-keep]
 # testtree <- drop.tip(phylo_a,phylo_a$tip.label[drop])
 # rm(drop, keep)

trees <- c("phylo_a_a_a", "phylo_b_a_a") # c("testtree", "phylo_a", "phylo_a_a_a", "phylo_a_b_a", "phylo_b", "phylo_b_a_a", "phylo_b_b_a")

for(i in 1:length(trees)){
  assign(paste("RD.", trees[i], sep=""), unlist(root.distance(get(trees[i]), mc.cores = 28)))
  #saveRDS(get(paste("RD.", trees[i], sep="")), paste("RD.", trees[i], ".rds", sep=""))
  print(paste(paste("RD.", trees[i], sep=""), "calculated.", Sys.time()))

  assign(paste("EDGES.", trees[i], sep=""), mclapply(1:Ntip(get(trees[i])), get_edges, phylo=get(trees[i]), mc.cores=4))
  #saveRDS(get(paste("EDGES.", trees[i], sep="")), paste("EDGES.", trees[i], ".rds", sep=""))
  print(paste(paste("EDGES.", trees[i], sep=""), "calculated.", Sys.time()))
}
rm(i)

####################
# 2. Calculate MRD #
####################

analyses <- cbind(
  c("phylo_a_a_a", "phylo_b_a_a"),
  c("MATCHES_a_a_a", "MATCHES_b_a_a")
)
rownames(analyses) <- c("a_a_a", "b_a_a")

# stub function for parallelization
mrd <- function(i, MATCHES, phylo, RD){
  MRD <- mean(RD[which(phylo$tip.label %in% MATCHES[MATCHES[,2] %in% colnames(comm)[comm[i,] == 1],1])])
  return(MRD)
}

for(n in rownames(analyses)){

  # choose relevant matching
  MATCHES <- get(analyses[n,2])
  phylo   <- get(analyses[n,1])
  RD      <- get(paste("RD.", analyses[n,1], sep=""))
  print(paste("Starting calculation of analysis", n))

  # calculate observed MRD
  MRD <- unlist(mclapply(1:nrow(comm), mrd, MATCHES=MATCHES, phylo=phylo, RD=RD, mc.cores = 28))
  print("observed MRD calculated")

  # calculate randomized MRDs
  MATCHES.rnd <- MATCHES
  for(i in 1:99){
    MATCHES.rnd[!is.na(MATCHES.rnd[,2]),2] <- sample(MATCHES.rnd[!is.na(MATCHES.rnd[,2]),2])
    MRD <- cbind(MRD, unlist(mclapply(1:nrow(comm), mrd, MATCHES=MATCHES.rnd, phylo=phylo, RD=RD, mc.cores = 28)))
    print(paste("replicate", i, "finished"))
  }

  colnames(MRD) <- c("obs", paste("rnd.", 1:99, sep=""))
  rownames(MRD) <- rownames(comm)

  assign(paste("MRD.", n, sep=""), as.data.frame(MRD))
}
rm(mrd, MRD, phylo, MATCHES, MATCHES.rnd, RD, i, n)

#######################
# 2. Compute phylosim #
#######################

for(n in rownames(analyses)){
  #n = rownames(analyses)[x]

  MATCHES <- get(analyses[n,2])
  phylo   <- get(analyses[n,1])
  EDGES   <- get(paste("EDGES.", analyses[n,1], sep=""))
  print(paste("Starting calculation of analysis", n))
  
  # reduce community matrix to species included in the relevant matching
  comm_red <- comm[,colnames(comm) %in% MATCHES[,2]]
  
  # rename community matrix with phylogeny tip labels
  comm.obs <- comm_red
  MATCHES.tmp <- MATCHES[!is.na(MATCHES[,2]),]
  idx <- as.vector(MATCHES.tmp$tip)
  names(idx) <- as.vector(MATCHES.tmp[,2])
  rm(MATCHES.tmp)
  colnames(comm.obs) <- idx[colnames(comm.obs)]
  rm(idx)
  
  ps_no <- list()
  ps_bl <- list()
  ps <- phylosim(phylo, comm.obs, EDGES = EDGES)
  ps_no[[1]] <- ps[[1]]
  ps_bl[[1]] <- ps[[2]]
  rm(ps)
  print("observed ps calculated")
 
  # calculate randomised phylosim
  
  for(i in 1:99){
    
    # shuffle matching
    MATCHES.rnd <- MATCHES
    MATCHES.rnd[!is.na(MATCHES.rnd[,2]),2] <- sample(MATCHES.rnd[!is.na(MATCHES.rnd[,2]),2])
    MATCHES.tmp <- MATCHES.rnd[!is.na(MATCHES.rnd[,2]),]
    idx <- as.vector(MATCHES.tmp$tip)
    names(idx) <- as.vector(MATCHES.tmp[,2])
    rm(MATCHES.tmp)
    # create new (randomized) community matrix
    comm.rnd <- comm_red
    colnames(comm.rnd) <- idx[colnames(comm.rnd)]
    rm(idx)
    
    ps <- phylosim(phylo, comm.rnd, EDGES = EDGES)
    ps_no[[i+1]] <- ps[[1]]
    ps_bl[[i+1]] <- ps[[2]]
    print(paste("replicate", i, "finished"))
    rm(comm.rnd, ps)
  }
  
  rm(i, comm_red, comm.obs, MATCHES, phylo, EDGES)
  assign(paste("ps_no.", n, sep=""), ps_no)
  assign(paste("ps_bl.", n, sep=""), ps_bl)
  #saveRDS(get(paste("ps.", n, sep="")), paste("ps.", n, ".rds", sep=""))
  rm(ps_no, ps_bl)
}
rm(n)

# for(n in rownames(analyses)){
#   assign(paste("ps.", n, sep=""), readRDS(paste("ps.", n, ".rds", sep="")))
# }
# rm(n)

save.image("phylostruct.RData")

save(comm, MRD.a_a_a, MRD.b_a_a, ps_no.a_a_a, ps_no.b_a_a, ps_bl.a_a_a, ps_bl.b_a_a, file="MRD_ps.RData")
