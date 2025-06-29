library(nnet)
real_dat <- readRDS("data/clean_data.RData")
K <- length(unique(real_dat$arm))

prpn_mod <- multinom(arm ~ baseline + gender + severity, 
                     subset(real_dat, select = c(arm, baseline, gender,
                                                 severity)))
prpn_est <- as.data.frame(predict(prpn_mod, type = "probs"))

twins_ind <- matrix(NA, nrow(real_dat), K)
colnames(twins_ind) <- paste("Arm", 1:K)
twins_ind <- as.data.frame(twins_ind)
twins_eff <- twins_ind
twins_safe <- twins_ind

for(i in 1:nrow(real_dat)){
  assgn_arm <- as.character(real_dat$arm[i])
  for(arm in paste("Arm", 1:K)){
    if(arm != assgn_arm){
      eligible_twins_ind <- which(real_dat$arm == arm)
      diff <- abs(prpn_est[[arm]][eligible_twins_ind] - prpn_est[[arm]][i])
      
      if(length(eligible_twins_ind) > 0){
        ind <- eligible_twins_ind[which.min(diff)]
        twins_ind[i, arm] <- ind
        twins_eff[i, arm] <- real_dat[["pct_change_salt"]][ind]
        twins_safe[i, arm] <- real_dat[["pct_change_lymph"]][ind]
      }
    } else{
      twins_ind[i, arm] <- i
      twins_eff[i, arm] <- real_dat[["pct_change_salt"]][i]
      twins_safe[i, arm] <- real_dat[["pct_change_lymph"]][i]
    }
  }
}

rm(list = ls())