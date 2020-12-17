##### Wrapper to extract level effects #####
# fit = model fitted
# level = which level, numeric
# trait = trait name
# ancestral = ancestral range
extractLEVELS_STAN <- function(fit, level, trait, ancestral) {
  efts <- ranef(fit) 
  efts_lvl <- efts[[level]]
  Intercepts <- data.frame(efts_lvl[ , , "Intercept"])
  Intercepts$Variable <- "Intercept"
  Intercepts$Realm <- rownames(Intercepts)
  Bio1 <- data.frame(efts_lvl[ , , "CHELSA_bio10_1"])
  Bio1$Variable <- "MAT"
  Bio1$Realm <- rownames(Bio1)
  Bio12 <- data.frame(efts_lvl[ , , "CHELSA_bio10_12"])
  Bio12$Variable <- "MAP"
  Bio12$Realm <- rownames(Bio12)
  Bio4 <- data.frame(efts_lvl[ , , "CHELSA_bio10_4"])
  Bio4$Variable <- "T Season"
  Bio4$Realm <- rownames(Bio4)
  Bio15 <- data.frame(efts_lvl[ , , "CHELSA_bio10_15"])
  Bio15$Variable <- "P Season"
  Bio15$Realm <- rownames(Bio15)
  RAD <- data.frame(efts_lvl[ , , "rad"])
  RAD$Variable <- "Rad" 
  RAD$Realm <- rownames(RAD)
  pH <- data.frame(efts_lvl[ , , "pH"])
  pH$Variable <- "pH"
  pH$Realm <- rownames(pH)
  CLY <- data.frame(efts_lvl[ , , "cly"])
  CLY$Variable <- "Cly"
  CLY$Realm <- rownames(CLY)
  DR <- data.frame(efts_lvl[ , , "DR"])
  DR$Variable <- "DR"
  DR$Realm <- rownames(DR)
  RD <- data.frame(efts_lvl[ , , "MRD"])
  RD$Variable <- "RD" 
  RD$Realm <- rownames(RD)
  
  levelEffects <- rbind(Intercepts, Bio1, Bio12, Bio4, Bio15, 
                        RAD, pH, CLY, DR, RD)
  levelEffects$Trait <- trait
  levelEffects$Ancestral <- ancestral
  return(levelEffects)
  
}

##### Wrapper that allow creating a final data.frame from the samples #####
#pars = "^[sd,b]_*"
#Traits <- c("SLA", "LNM", "LPM")
#Ancestral <- c("Arid", "Equatorial", "Temperate", "Polar", "Tundra", "Other")

extractSAMPLES_Stan <- function(fit, pars, ancestral, trait) {
  sampOBJ <- posterior_samples(fit, pars = pars)
  sampOBJ$Trait <- trait
  sampOBJ$Ancestral <- ancestral
  return(sampOBJ)
}

##### Wrapper to extract fixed effects #####
# fit = model fitted
# trait = trait name
# ancestral = ancestral range
extractFIXED_STAN <- function(fit, trait, ancestral) {
  efts <- data.frame(fixef(fit)) 
  efts$Covars <- rownames(efts)
  efts$Trait <- trait
  efts$Ancestral <- ancestral
  
  return(efts) 
}

##### Wrapper to extract fitted values #####

extractFITTED_STAN <- function(fit, probs, trait, ancestral) {
  fitVal <- data.frame(fitted(equa_lpm_fit, summary = TRUE, 
                              robust = TRUE, probs = probs))
  fitVal$Trait <- trait
  fitVal$Ancestral <- ancestral
  
  return(fitVal)
}

##### Wrapper to extract and combine fitted values with data #####

extract_FITTED_DATA_STAN <-  function(fit, probs, trait, ancestral) {
  fitVal <- data.frame(fitted(fit, summary = TRUE, 
                              robust = TRUE, probs = probs))
  fitVal$Trait <- trait
  fitVal$Ancestral <- ancestral
  DATA <- fit$data
  
  fitDATA <- cbind(fitVal, DATA)
  
  return(fitDATA)

}

##### Wrapper to extract and combine Predicted values with data #####

extract_PREDICTED_DATA_STAN <-  function(fit, probs, trait, ancestral) {
  predVal <- data.frame(predict(fit, summary = TRUE, 
                              robust = TRUE, probs = probs))
  predVal$Trait <- trait
  predVal$Ancestral <- ancestral
  DATA <- fit$data
  
  predDATA <- cbind(predVal, DATA)
  
  return(predDATA)
  
}
