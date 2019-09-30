# first argument: tree to be used ("a" or "b")
# second argument: which matching to use ("a" or "b")
args = commandArgs(trailingOnly=TRUE)
if(length(args)!=2) stop("Please provide the correct arguments.")

tree = args[1]
matching = args[2]

#setwd("~/Documents/WOLF/PROJECTS/58 World Checklist paper/analyses 2019")
setwd("/data_vol/wolf/WCSP_paper")

library(ape)
library(phytools)

load("MATCHES.RData")

source("functions.R")

########################
# Analysis starts here #
########################

# get list of all accepted speices
goodspp <- wcsp[wcsp$genus_hybrid_marker == "" & wcsp$species_hybrid_marker == "" & wcsp$infraspecific_rank == "" & wcsp$species != "" & wcsp$taxon_status_description == "Accepted",]

# Tree
phylo <- get(paste("phylo_", tree, sep=""))
# Matching
MATCHES <- get(paste("MATCHES_", tree, "_", matching, sep=""))

# get list of good species that are not in matching
toadd <- goodspp[!goodspp$checklist_id %in% MATCHES[!is.na(MATCHES[,2]),2],]

# initialize vector to record which method to use to find MRCA (genus = 1, family = 2, or order = 3)
method <- vector("numeric", nrow(toadd))

# get a list of genera that are already represented in the tree
genera <- unique(as.vector(goodspp[goodspp$checklist_id %in% as.vector(MATCHES[,2]),"genus"]))

# species that can be added to genus
method[toadd$genus %in% genera] <- 1

# get the families of the species that remain to be added
families <- as.vector(unique(toadd$accepted_family))

# species that can be added to family
method[method==0][as.vector(toadd[method==0,"accepted_family"]) %in% as.vector(goodspp[goodspp$checklist_id %in% as.vector(MATCHES[,2]),"accepted_family"])] <- 2

# species that can be added to order
method[method==0][!as.vector(toadd[method==0,"accepted_family"]) %in% as.vector(goodspp[goodspp$checklist_id %in% as.vector(MATCHES[,2]),"accepted_family"])] <- 3

# create index for orders
{unresolved_families <- list()
unresolved_families[["Cynomoriaceae"]] <- c("Altingiaceae", "Aphanopetalaceae", "Cercidiphyllaceae", "Crassulaceae", "Daphniphyllaceae", "Grossulariaceae", "Haloragaceae", "Hamamelidaceae", "Iteaceae", "Paeoniaceae", "Penthoraceae", "Peridiscaceae", "Saxifragaceae", "Tetracarpaeaceae")
unresolved_families[["Circaeasteraceae"]] <- c("Berberidaceae", "Eupteleaceae", "Lardizabalaceae", "Menispermaceae", "Papaveraceae", "Ranunculaceae")
unresolved_families[["Mitrastemonaceae"]] <- c("Actinidiaceae", "Balsaminaceae", "Cyrillaceae", "Clethraceae", "Diapensiaceae", "Ebenaceae", "Ericaceae", "Fouquieriaceae", "Lecythidaceae", "Marcgraviaceae", "Pentaphylacaceae", "Polemoniaceae", "Primulaceae", "Roridulaceae", "Sapotaceae", "Sarraceniaceae", "Sladeniaceae", "Styracaceae", "Symplocaceae", "Tetrameristaceae", "Theaceae")
unresolved_families[["Cytinaceae"]] <- c("Bixaceae", "Cistaceae", "Dipterocarpaceae", "Malvaceae", "Muntingiaceae", "Neuradaceae", "Sarcolaenaceae", "Sphaerosepalaceae", "Thymelaeaceae")
unresolved_families[["Physenaceae"]] <- c("Achatocarpaceae", "Aizoaceae", "Amaranthaceae", "Anacampserotaceae", "Ancistrocladaceae", "Asteropeiaceae", "Barbeuiaceae", "Basellaceae", "Cactaceae", "Caryophyllaceae", "Didiereaceae", "Dioncophyllaceae", "Droseraceae", "Drosophyllaceae", "Frankeniaceae", "Gisekiaceae", "Halophytaceae", "Kewaceae", "Limeaceae", "Lophiocarpaceae", "Macarthuriaceae", "Microteaceae", "Molluginaceae", "Montiaceae", "Nepenthaceae", "Nyctaginaceae", "Phytolaccaceae", "Plumbaginaceae", "Polygonaceae", "Portulacaceae", "Rhabdodendraceae", "Petiveriaceae", "Sarcobataceae", "Simmondsiaceae", "Stegnospermataceae", "Talinaceae", "Tamaricaceae")
unresolved_families[["Tetracarpaeaceae"]] <- c("Altingiaceae", "Aphanopetalaceae", "Cercidiphyllaceae", "Crassulaceae", "Cynomoriaceae", "Daphniphyllaceae", "Grossulariaceae", "Haloragaceae", "Hamamelidaceae", "Iteaceae", "Paeoniaceae", "Penthoraceae", "Peridiscaceae", "Saxifragaceae")}

# reduce toadd for testing purposes
# toadd <- toadd[1:100,]

for(i in 1:nrow(toadd)){
  
  # define new taxon name
  newname <- gsub(pattern = " ", replacement = "_", toadd[i,"accepted_name"])
  # if name already exists in tree, append "_WLE" to make it unique
  if(newname %in% phylo$tip.label) newname <- paste(newname, "_WLE", sep="")
  
  # gather the tips representing the higher taxon of interest
  if(method[i] == 1){ # from genus
    tips = MATCHES[MATCHES[,2] %in% as.vector(goodspp[goodspp$genus == as.vector(toadd[i,"genus"]),"checklist_id"]),1]
  }
  if(method[i] == 2){ # from family
    tips = MATCHES[MATCHES[,2] %in% as.vector(goodspp[goodspp$accepted_family == as.vector(toadd[i,"accepted_family"]),"checklist_id"]),1]
  }
  if(method[i] == 3){ # from order
    tips = MATCHES[MATCHES[,2] %in% as.vector(goodspp[goodspp$accepted_family %in% unresolved_families[[as.vector(toadd[i,"accepted_family"])]],"checklist_id"]),1]
  }
  
  if(length(tips) < 1) stop("Attempting to add to a group that is not represented in the tree!")
  
  if(length(tips) == 1){ # if adding to a terminal branch
    nn <- which(phylo$tip.label == tips)
    phylo <- bind.tip2(phylo, newname, where = nn, position = 0.5 * phylo$edge.length[which(phylo$edge[, 2] == nn)])
  } else { # if adding to an internal node
    nn <- findMRCA(phylo, tips = tips)
    phylo <- bind.tip2(phylo, newname, where = nn)
  }
  MATCHES <- rbind(MATCHES, c(newname, as.vector(toadd[i,"checklist_id"])))

  if(i %% ceiling(nrow(toadd)/100) == 0) print(paste(i/ceiling(nrow(toadd)/100), "% complete", Sys.time()))
}

phylo <- untangle(phylo)

assign(paste("phylo_", tree, "_", matching, "_a", sep=""), phylo)
assign(paste("MATCHES_", tree, "_", matching, "_a", sep=""), MATCHES)

save(list=c(paste("phylo_", tree, sep=""), paste("phylo_", tree, "_", matching, "_a", sep=""), paste("MATCHES_", tree, "_", matching, sep=""), paste("MATCHES_", tree, "_", matching, "_a", sep=""), "wcsp"), file = paste("MATCHES_", tree, "_", matching, ".RData", sep=""))

# keep <- sample(1:length(phylo$tip.label), 300)
# drop <- 1:length(phylo$tip.label)
# drop <- drop[-keep]
# testtree <- drop.tip(phylo,phylo$tip.label[drop])
# 
# testmatch <- MATCHES[keep,]
# rownames(testmatch) <- testmatch$tip
# testmatch <- testmatch[testtree$tip.label,]  
# 
# write.tree(testtree, "testttee.tre")
# 
# nn <- findMRCA(testtree, tips = testtree$tip.label[testmatch$conservative %in% as.vector(goodspp[goodspp$accepted_family == as.vector(add_to_family_multi[i,"accepted_family"]),"checklist_id"])])
# 
# testtree2 <- bind.tip(testtree, "test", where = nn)
# 
# write.tree(testtree2, "testttee2.tre")