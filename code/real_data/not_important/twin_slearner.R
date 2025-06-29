library(randomForest)
real_dat <- readRDS("data/clean_data.RData")
K <- length(unique(real_dat$arm))
rf_eff_mods <- vector(mode = "list", length = K)
rf_safe_mods <- rf_eff_mods
countfact_eff <- matrix(NA, nrow(real_dat), K)
colnames(countfact_eff) <- paste("Arm", 1:K)
countfact_eff <- as.data.frame(countfact_eff)
countfact_safe <- countfact_eff

for(k in 1:K){
  assgn_arm <- paste("Arm", k)
  rf_eff_mods[[k]] <- randomForest(pct_change_salt ~ ., 
                               subset(real_dat, subset = arm == assgn_arm, 
                                      select = -c(arm, pct_change_lymph))
                               #, ntree = 100, maxnodes = 6
                               )
  rf_safe_mods[[k]] <- randomForest(pct_change_lymph ~ ., 
                                   subset(real_dat, subset = arm == assgn_arm, 
                                          select = -c(arm, pct_change_salt))
                                   #, ntree = 100, maxnodes = 6
                                    )
}

for(k in 1:K){
  assgn_arm <- paste("Arm", k)
  pred_eff <- predict(rf_eff_mods[[k]], 
                      subset(real_dat, subset = arm != assgn_arm,
                             select = -c(arm, pct_change_lymph, 
                                         pct_change_salt)))
  countfact_eff[real_dat$arm != assgn_arm, assgn_arm] <- pred_eff
  countfact_eff[real_dat$arm == assgn_arm, assgn_arm] <- 
    real_dat[["pct_change_salt"]][real_dat$arm == assgn_arm]
  
  pred_safe <- predict(rf_safe_mods[[k]], 
                      subset(real_dat, subset = arm != assgn_arm,
                             select = -c(arm, pct_change_lymph, 
                                         pct_change_salt)))
  countfact_safe[real_dat$arm != assgn_arm, assgn_arm] <- pred_safe
  countfact_safe[real_dat$arm == assgn_arm, assgn_arm] <- 
    real_dat[["pct_change_lymph"]][real_dat$arm == assgn_arm]
  
}

saveRDS(countfact_eff, "output/real_dat/countfact_eff.RData")
saveRDS(countfact_safe, "output/real_dat/countfact_safe.RData")

rm(list = ls())