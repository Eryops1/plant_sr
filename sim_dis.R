#setwd("~/Documents/WOLF/PROJECTS/58 World Checklist paper/analyses 2019")

source("plant_sr/functions.R")

load("plant_sr_data/comm.RData")

  bs_spp <- betasim(comm)
  
  save(bs_spp, file = "bs_spp.RData")
  
  save.image("sim_dis.RData")
