setwd("C:/Users/jpintole/Dropbox/Lucie_Jesus/RCN_evolutionary_legacy")

library(ape)
library(phytools)
library(dplyr)
library(tidyr)
library(picante)
library(brms)
library(cmdstanr)

set_cmdstan_path(path = NULL)

halfScale <- function(x){
  val <- (x-mean(x)) * 0.5/sd(x)
  return(val)
}

##### Load data #####
traits <- read.csv("Data/TryData/FinalDATA/traitsMATCH_checkNOnas_FINAL.csv")
names(traits)

phylo <- multi2di(read.nexus("Data/Phylogeny/phylo_EL.nex"))
is.binary.phylo(phylo)
is.ultrametric(phylo)

phyloMatch <- drop.tip(phylo, setdiff(phylo$tip.label, traits$Species))

##### Subsets by ancestral range #####
### Ancestral ranges based on DEC-J model
# A = Equatorial
# B = Arid
# C = Temperate
# D = Snow/Tundra
# E = Polar

traits %>% 
  group_by(Major_climate_KG_family) %>% 
  count()

equatorial <- droplevels(subset(traits, Major_climate_KG_family == "A"))
equaPhy <- drop.tip(phyloMatch, setdiff(phyloMatch$tip.label, equatorial$Species))

##### Run analysis by functional trait and by ancestral range #####
dir.create("Multilevel/ANCESTRAL_MODELS/Equatorial")

## Scale environmental covariates 
head(equatorial[, 20:26])

equatorial[, 20:26] <- apply(equatorial[, 20:26], 2, FUN = scale)
head(equatorial[, 20:26])

equatorial$DR <- log(equatorial$DR+1)
equatorial$MRD <- log(equatorial$MRD+1)

## Variance-covariance matrix
A_EQUA <- ape::vcv.phylo(equaPhy, corr = TRUE)

equatorial %>% 
  count(Realm, Site)

##### MODELS with interactions between covariates #####
traits <- c("SLA", "LNM", "LPM")

names(equatorial)

### Inits
level <- names(equatorial[c(3, 4, 11)])
covars <- names(equatorial[c(20:26, 37, 39)])
traitsMODEL <- names(equatorial[c(12, 13, 14)])

### Run Bayesian Phylogenetic Multilevel Model 
# This take a lot of time to run, please be patient!

for(i in 1:length(traits)){
  print(traits[i])
  DATA <- cbind(equatorial[, traitsMODEL[i]], equatorial[, covars], equatorial[, level])
  names(DATA) <- c("trait", names(DATA[2:13]))

  if (!file.exists(paste0("Multilevel/ANCESTRAL_MODELS/Equatorial/", traits[i], 
                          "_EQUATORIAL_4Chains.rds", sep = ""))) {
    equa_MODEL <- brm(
      log(trait) ~ 
        CHELSA_bio10_1 + CHELSA_bio10_12 + CHELSA_bio10_4 + CHELSA_bio10_15 + 
        rad + pH + cly +
        DR + MRD + 
        DR*CHELSA_bio10_1 + DR*CHELSA_bio10_12 + DR*CHELSA_bio10_4 + DR*CHELSA_bio10_15 + 
        DR*rad + DR*pH + DR*cly +  
        MRD*CHELSA_bio10_1 + MRD*CHELSA_bio10_12 + MRD*CHELSA_bio10_4 + MRD*CHELSA_bio10_15 + 
        MRD*rad + MRD*pH + MRD*cly + 
        (1 | Phylogeny) + (1 | Species) +  
        (1 + CHELSA_bio10_1 + CHELSA_bio10_12 + CHELSA_bio10_4 + CHELSA_bio10_15 + 
           rad + pH + cly + DR + MRD | Realm),
      data = DATA,
      cov_ranef = list(Phylogeny = A_EQUA),
      chains = 4, iter = 10000, warmup = 2000, cores = getOption("mc.cores", 32), silent = TRUE,
      control = list(max_treedepth = 12, adapt_delta = 0.9999)
    )
    
    saveRDS(equa_MODEL, paste0("Multilevel/ANCESTRAL_MODELS/Equatorial/", traits[i], 
                               "_EQUATORIAL_4Chains_rep.rds", sep = ""))
  }
  
  print(paste0("Congratulations your Bayesian Phylogenetic Multilevel Model for ", traits[i], 
               " is complete, please check folder ANCESTRAL_MODELS/Equatorial", sep = ""))
}

