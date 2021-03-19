
##### Inits

#fites <- c("LPM_TUNDRA_POLAR_4Chains", "LPM_TUNDRA_noPhylo_4Chains")
#address <- "Dropbox/Lucie_Jesus/RCN_evolutionary_legacy/Multilevel/ANCESTRAL_MODELS/Tundra-Polar/"
#N = 2

# Wrapper to run Bayesian model selection using LOO - Efficient approximate leave-one-out cross-validation
extractLOO <- function(fits, direct, nModels) { 
 
   if ( ! ("loo" %in% installed.packages())) {install.packages("loo", dependencies = T)} 
  
  require(loo)
  loos <- list()
  
  for(i in 1:nModels) { 
    tm <- fits[i]
    
    fit <- readRDS(paste0(direct, tm, ".rds"))
    
    loos[[i]] <- loo(fit)
    
    print(paste0("loo estimated for model ", tm))
  }
  
  comparison <- loo_compare(x = loos)
  return(comparison)
}

## Usage
# The first row correspond to the best model
#ss <- extractLOO(fits = fites, direct = address, nModels = N)
#print(ss, simplify = FALSE, digits = 3) 
