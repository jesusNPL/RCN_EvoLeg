setwd("Dropbox/Lucie_Jesus/RCN_evolutionary_legacy")

library(ape)
library(phytools)
library(dplyr)
library(tidyr)
library(picante)

traits <- read.csv("Data/TryData/gapFilling/traitsMATCH_checkNOnas.csv")
names(traits)

sppTraits <- as.character(unique(traits$correct.names))

phylo <- multi2di(read.nexus("Data/Phylogeny/phylo_EL.nex"))
is.binary.phylo(phylo)
is.ultrametric(phylo)

phyloMatch <- drop.tip(phylo, setdiff(phylo$tip.label, sppTraits))

sppPhylo <- phyloMatch$tip.label

traitsMatch <- traits[traits$correct.names %in% sppPhylo, ]
length(unique(traitsMatch$correct.names))

traitsMatch$Species_model <- traitsMatch$correct.names # species efect
traitsMatch$Species <- traitsMatch$correct.names # phylogenetic effect

rm(phylo, traits)

traitsMatch %>% group_by(Major_climate_KG_family) %>% summarise(non_na_count = sum(!is.na(sla)))
traitsMatch %>% group_by(Major_climate_KG_family) %>% summarise(non_na_count = sum(!is.na(lnm)))
traitsMatch %>% group_by(Major_climate_KG_family) %>% summarise(non_na_count = sum(!is.na(lpm)))

traitsMatch %>% group_by(Major_climate_KG_family) %>% summarise(non_na_count = sum(!is.na(sla_imp)))
traitsMatch %>% group_by(Major_climate_KG_family) %>% summarise(non_na_count = sum(!is.na(lnm_imp)))
traitsMatch %>% group_by(Major_climate_KG_family) %>% summarise(non_na_count = sum(!is.na(lpm_imp)))

traitsMatch <- traitsMatch[, c("Species", "Species_model", "genus", "family", "PFTs",  
                               "ObservationID", "LON_site", "LAT_site", "Realm", # taxonomic and site variables
                               "sla", "lnm", "lpm", "lnp",  # Traits original
                               "sla_imp", "lnm_imp", "lpm_imp", "lnp_imp", # Traits imputed 
                               "CHELSA_bio10_1", "CHELSA_bio10_4", "CHELSA_bio10_12", "CHELSA_bio10_15", # Climatic covariates
                               "rad", "pH", "cly", # Physical and soil covariates  
                               "Ancestral_DECJ_family", "Major_climate_KG_family", # Ancestral range information
                               "KG_classes", "KG_class_broad", "KG_class_broader")]
traitsMatch$Major_climate_KG_family <- factor(traitsMatch$Major_climate_KG_family, 
                                              levels = c("E", "D", "C", "B", "A"))
traitsMatch %>% 
  count(Realm, Major_climate_KG_family)

#traitsMatch2 <- traitsMatch %>% tidyr::drop_na(Realm)
Sites <- traitsMatch %>% 
  group_by(LON_site, LAT_site) %>%  
  summarize(longitude = mean(LON_site), 
            latitude = mean(LAT_site))

coms <- rep("site", nrow(Sites))
lst <- NULL

for(i in 1:nrow(Sites)) {
  lst[i] <- paste0(coms[i], "_", i)
}

Sites$Site <- lst
head(Sites)

traitsMatch <- full_join(traitsMatch, Sites, by = c("LON_site", "LAT_site"))

SR <- traitsMatch %>% 
  group_by(Site) %>% 
  count(Species) %>% 
  summarise(SR = sum(n))

tail(SR)

traitsMatch <- full_join(traitsMatch, SR , by = "Site")

unique(traitsMatch$Site)


comus <- traitsMatch %>% 
  count(Site) %>% 
  filter(n >= 2)

x <- SR %>% 
  filter(SR >= 2)

##### Diversification estimation #####
eq <- evol.distinct(phyloMatch, type = "equal.splits", use.branch.lengths = TRUE)

eq$DR <- 1/eq$w
head(eq)

##### Get species ages #####
## first get the node numbers of the tips
nodes <- sapply(phyloMatch$tip.label, function(x, y) which(y == x), y = phyloMatch$tip.label)
## then get the edge lengths for those nodes
edge.lengths <- setNames(phyloMatch$edge.length[sapply(nodes, function(x, y) which(y==x), 
                                                       y = phyloMatch$edge[, 2])], names(nodes))
ageSPP <- as.data.frame(edge.lengths)
ageSPP$Species <- rownames(ageSPP)
names(ageSPP) <- c("spp_AGE", "Species")
head(ageSPP)

##### Get mean root distance #####
tree <- ape::compute.brlen(phyloMatch, 1)
allDists <- ape::dist.nodes(tree)

#subset allDists to those between the root and the tips
rootDists <- allDists[length(tree$tip.label) + 1, 1:length(tree$tip.label)]

#bind these root distances into a dataframe with species labels
tipsToRoot <- data.frame(Species = tree$tip.label, Nnodes = rootDists)
head(tipsToRoot)

### Join species age, mean root distance and diversification rate 
DR_AGE <- dplyr::full_join(eq, ageSPP, by = "Species")
names(DR_AGE) <- c("Species", "EqualSplit", "DR", "Age")
head(DR_AGE)

DR_AGE_MRD <- dplyr::full_join(DR_AGE, tipsToRoot, by = "Species")
names(DR_AGE_MRD) <- c("Species", "EqS", "DR", "Age", "MRD")
head(DR_AGE_MRD)

cor.test(DR_AGE_MRD$DR, DR_AGE_MRD$Age)
cor.test(DR_AGE_MRD$DR, DR_AGE_MRD$MRD)
cor.test(DR_AGE_MRD$Age, DR_AGE_MRD$MRD)

##### FINAL JOIN #####

traitsMatch <- full_join(traitsMatch, DR_AGE_MRD, by = "Species")

cor.test(traitsMatch$SR, traitsMatch$MRD)

write.csv(traitsMatch, file = "Data/TryData/FinalDATA/traitsMATCH_checkNOnas_FINAL.csv")
