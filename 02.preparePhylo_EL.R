

library(V.PhyloMaker)
library(ape)

source("https://raw.githubusercontent.com/jesusNPL/ManageTRY/master/demonCheckScinames.R")

dato <- read.csv("Data/TryData/traitsMATCH_checkNOnas_FINAL.csv")

scinames <- sort(as.character(unique(dato$correct.names)))
scinames <- gsub("_", " ", scinames)

###### Check scientific names #####
scinamesTPL <- check_TPLScinames(scinames = scinames, saveResults = FALSE)

write.csv(scinamesTPL, "Data/Phylogeny/traits_TPL_check.csv")

##### Prepare phylogeny using V.Phylomaker
data_traits <- data.frame(species = scinamesTPL$TPLSciname, 
                                 genus = scinamesTPL$TPLGenus, 
                                 family = scinamesTPL$TPLFamily)

phylo <- phylo.maker(data_traits, scenarios = "S3", output.tree = TRUE)

write.nexus(phylo$scenario.3, file = "Data/Phylogeny/phylo_EL.nex")
write.nexus(phylo$tree.scenario.3, file = "Data/Phylogeny/phylo_BackBone.nex")
write.csv(phylo$species.list, "Data/Phylogeny/data_used_phylogeny.csv")
