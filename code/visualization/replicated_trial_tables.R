file_choice <- 3
ts_sim <- readRDS(paste("output/ts_sim_", file_choice, ".RData", sep = ""))
rand_sim <- readRDS(paste("output/rand_sim_", 1, ".RData", sep = ""))
rits_sim <- readRDS(paste("output/rits_sim_", file_choice, ".RData", sep = ""))
n_iter <- length(ts_sim)
N <- length(ts_sim[[1]]$trt)
ind <- c(seq(50, 100, 10), 150, 200)
K <- length(unique(rand_sim[[1]]$trt))

# # average cumulative regret (both continuous end-points)
# tab_metric <- matrix(0, 6, length(ind))
# 
# tmp <- t(sapply(1:n_iter, function(i){
#   cumsum(rand_sim[[i]]$regret)[ind]
# }))
# tab_metric[1, ] <- apply(tmp, 2, mean)
# 
# tmp <- t(sapply(1:n_iter, function(i){
#   cumsum(ts_sim[[i]]$regret)[ind]
# }))
# tab_metric[2, ] <- apply(tmp, 2, mean)
# 
# tmp <- t(sapply(1:n_iter, function(i){
#   cumsum(rits_sim[[i]]$regret)[ind]
# }))
# tab_metric[3, ] <- apply(tmp, 2, mean)
# 
# tmp <- t(sapply(1:n_iter, function(i){
#   cumsum(rand_sim[[i]]$subopt_trt)[ind]
# }))
# tab_metric[4, ] <- apply(tmp, 2, mean)
# 
# tmp <- t(sapply(1:n_iter, function(i){
#   cumsum(ts_sim[[i]]$subopt_trt)[ind]
# }))
# tab_metric[5, ] <- apply(tmp, 2, mean)
# 
# tmp <- t(sapply(1:n_iter, function(i){
#   cumsum(rits_sim[[i]]$subopt_trt)[ind]
# }))
# tab_metric[6, ] <- apply(tmp, 2, mean)

### Quality of Causal Effects

sim_dat <- readRDS("metadata/sim_dat.RData")
mu_true <- sim_dat$mu_true * (1/5)
# mu_true <- apply(X %*% beta_true[, , 1], 2, mean)
contr_true <- mu_true - mu_true[1]
contr_true <- contr_true[-1]

sim_choice <- readRDS("metadata/sim_choice.RData")
ate_start <- sim_choice$ate_start
ate_ind <- ind - ate_start
std_ind <- ind - 49

# RAND
sim <- rand_sim
covr_ate_rand <- matrix(0, K, length(ind))
width_ate_rand <- matrix(0, K, length(ind))
bias_ate_rand <- matrix(0, K, length(ind))
covr_contr_rand <- matrix(0, K-1, length(ind))
width_contr_rand <- matrix(0, K-1, length(ind))
bias_contr_rand <- matrix(0, K-1, length(ind))
rmse_contr_rand <- matrix(0, K-1, length(ind))
covr_contr_std <- matrix(0, K-1, length(ind))
width_contr_std <- matrix(0, K-1, length(ind))
bias_contr_std <- matrix(0, K-1, length(ind))
rmse_contr_std <- matrix(0, K-1, length(ind))
# stop_n_rand <- sapply(1:n_iter, function(i) stop_trial(sim[[i]], thres))
for(k in 1:K){
  for(i in 1:length(ind)){
    for(iter in 1:n_iter){
      cond <- (sim[[iter]]$ate[ate_ind[i], k, 2] <= mu_true[k]) &&
        (sim[[iter]]$ate[ate_ind[i], k, 3] >= mu_true[k])
      covr_ate_rand[k, i] <- covr_ate_rand[k, i] + ifelse(cond, 1, 0)
      width_ate_rand[k, i] <- width_ate_rand[k, i] + 
        abs(sim[[iter]]$ate[ate_ind[i], k, 3] - sim[[iter]]$ate[ate_ind[i], k, 2])
      bias_ate_rand[k, i] <- bias_ate_rand[k, i] + 
        sim[[iter]]$ate[ate_ind[i], k, 1] - mu_true[k]
      
      if(k > 1){
        cond <- (sim[[iter]]$contr[ate_ind[i], k-1, 2] <= contr_true[k-1]) &&
          (sim[[iter]]$contr[ate_ind[i], k-1, 3] >= contr_true[k-1])
        covr_contr_rand[k-1, i] <- covr_contr_rand[k-1, i] + ifelse(cond, 1, 0)
        width_contr_rand[k-1, i] <- width_contr_rand[k-1, i] + 
          abs(sim[[iter]]$contr[ate_ind[i], k-1, 3] - sim[[iter]]$contr[ate_ind[i], k-1, 2])
        bias_contr_rand[k-1, i] <- bias_contr_rand[k-1, i] + 
          sim[[iter]]$contr[ate_ind[i], k-1, 1] - contr_true[k-1]
        rmse_contr_rand[k-1, i] <- rmse_contr_rand[k-1, i] + 
          (sim[[iter]]$contr[ate_ind[i], k-1, 1] - contr_true[k-1])^2
        
        cond <- (sim[[iter]]$contr_standard[std_ind[i], k-1, 2] <= contr_true[k-1]) &&
          (sim[[iter]]$contr_standard[std_ind[i], k-1, 3] >= contr_true[k-1])
        covr_contr_std[k-1, i] <- covr_contr_std[k-1, i] + ifelse(cond, 1, 0)
        width_contr_std[k-1, i] <- width_contr_std[k-1, i] + 
          abs(sim[[iter]]$contr_standard[std_ind[i], k-1, 3] - 
                sim[[iter]]$contr_standard[std_ind[i], k-1, 2])
        bias_contr_std[k-1, i] <- bias_contr_std[k-1, i] + 
          sim[[iter]]$contr_standard[std_ind[i], k-1, 1] - contr_true[k-1]
        rmse_contr_std[k-1, i] <- rmse_contr_std[k-1, i] + 
          (sim[[iter]]$contr_standard[std_ind[i], k-1, 1] - contr_true[k-1])^2
      }
    }
    covr_ate_rand[k, i] <- covr_ate_rand[k, i] / n_iter
    width_ate_rand[k, i] <- width_ate_rand[k, i] / n_iter
    bias_ate_rand[k, i] <- bias_ate_rand[k, i] / n_iter
    
    if(k > 1){
      covr_contr_rand[k-1, i] <- covr_contr_rand[k-1, i] / n_iter
      width_contr_rand[k-1, i] <- width_contr_rand[k-1, i] / n_iter
      bias_contr_rand[k-1, i] <- bias_contr_rand[k-1, i] / n_iter
      rmse_contr_rand[k-1, i] <- sqrt(rmse_contr_rand[k-1, i] / (n_iter))
      
      covr_contr_std[k-1, i] <- covr_contr_std[k-1, i] / n_iter
      width_contr_std[k-1, i] <- width_contr_std[k-1, i] / n_iter
      bias_contr_std[k-1, i] <- bias_contr_std[k-1, i] / n_iter
      rmse_contr_std[k-1, i] <- sqrt(rmse_contr_std[k-1, i] / (n_iter))
    }
  }
}


# TS (Thompson Sampling with Benefit only)
sim <- ts_sim
covr_ate_ts <- matrix(0, K, length(ind))
width_ate_ts <- matrix(0, K, length(ind))
bias_ate_ts <- matrix(0, K, length(ind))
covr_contr_ts <- matrix(0, K-1, length(ind))
width_contr_ts <- matrix(0, K-1, length(ind))
bias_contr_ts <- matrix(0, K-1, length(ind))
rmse_contr_ts <- matrix(0, K-1, length(ind))
# stop_n_ts <- sapply(1:n_iter, function(i) stop_trial(sim[[i]], thres))
for(k in 1:K){
  for(i in 1:length(ind)){
    tmp <- which(sapply(1:n_iter, 
                        function(i) match(TRUE, sim[[i]]$ate[, k, 1] != 0)) <= ate_ind[i])
    for(iter in tmp){
      cond <- (sim[[iter]]$ate[ate_ind[i], k, 2] <= mu_true[k]) &&
        (sim[[iter]]$ate[ate_ind[i], k, 3] >= mu_true[k])
      covr_ate_ts[k, i] <- covr_ate_ts[k, i] + ifelse(cond, 1, 0)
      width_ate_ts[k, i] <- width_ate_ts[k, i] + 
        abs(sim[[iter]]$ate[ate_ind[i], k, 3] - sim[[iter]]$ate[ate_ind[i], k, 2])
      bias_ate_ts[k, i] <- bias_ate_ts[k, i] + 
        sim[[iter]]$ate[ate_ind[i], k, 1] - mu_true[k]
      
      if(k > 1){
        cond <- (sim[[iter]]$contr[ate_ind[i], k-1, 2] <= contr_true[k-1]) &&
          (sim[[iter]]$contr[ate_ind[i], k-1, 3] >= contr_true[k-1])
        covr_contr_ts[k-1, i] <- covr_contr_ts[k-1, i] + ifelse(cond, 1, 0)
        width_contr_ts[k-1, i] <- width_contr_ts[k-1, i] + 
          abs(sim[[iter]]$contr[ate_ind[i], k-1, 3] - sim[[iter]]$contr[ate_ind[i], k-1, 2])
        bias_contr_ts[k-1, i] <- bias_contr_ts[k-1, i] + 
          sim[[iter]]$contr[ate_ind[i], k-1, 1] - contr_true[k-1]
        rmse_contr_ts[k-1, i] <- rmse_contr_ts[k-1, i] + 
          (sim[[iter]]$contr[ate_ind[i], k-1, 1] - contr_true[k-1])^2
      }
    }
    covr_ate_ts[k, i] <- covr_ate_ts[k, i] / n_iter
    width_ate_ts[k, i] <- width_ate_ts[k, i] / n_iter
    bias_ate_ts[k, i] <- bias_ate_ts[k, i] / n_iter
    
    if(k > 1){
      covr_contr_ts[k-1, i] <- covr_contr_ts[k-1, i] / n_iter
      width_contr_ts[k-1, i] <- width_contr_ts[k-1, i] / n_iter
      bias_contr_ts[k-1, i] <- bias_contr_ts[k-1, i] / n_iter
      rmse_contr_ts[k-1, i] <- sqrt(rmse_contr_ts[k-1, i] / n_iter)
    }
  }
}


# RiTS (Thompson Sampling with Benefit and Binary Risk)
sim <- rits_sim
covr_ate_rits <- matrix(0, K, length(ind))
width_ate_rits <- matrix(0, K, length(ind))
bias_ate_rits <- matrix(0, K, length(ind))
covr_contr_rits <- matrix(0, K-1, length(ind))
width_contr_rits <- matrix(0, K-1, length(ind))
bias_contr_rits <- matrix(0, K-1, length(ind))
rmse_contr_rits <- matrix(0, K-1, length(ind))
# stop_n_rits <- sapply(1:n_iter, function(i) stop_trial(sim[[i]], thres))
for(k in 1:K){
  for(i in 1:length(ind)){
    tmp <- which(sapply(1:n_iter, 
                        function(i) match(TRUE, sim[[i]]$ate[, k, 1] != 0)) <= ate_ind[i])
    for(iter in tmp){
      cond <- (sim[[iter]]$ate[ate_ind[i], k, 2] <= mu_true[k]) &&
        (sim[[iter]]$ate[ate_ind[i], k, 3] >= mu_true[k])
      covr_ate_rits[k, i] <- covr_ate_rits[k, i] + ifelse(cond, 1, 0)
      width_ate_rits[k, i] <- width_ate_rits[k, i] + 
        abs(sim[[iter]]$ate[ate_ind[i], k, 3] - sim[[iter]]$ate[ate_ind[i], k, 2])
      bias_ate_rits[k, i] <- bias_ate_rits[k, i] + 
        sim[[iter]]$ate[ate_ind[i], k, 1] - mu_true[k]
      
      if(k > 1){
        cond <- (sim[[iter]]$contr[ate_ind[i], k-1, 2] <= contr_true[k-1]) &&
          (sim[[iter]]$contr[ate_ind[i], k-1, 3] >= contr_true[k-1])
        covr_contr_rits[k-1, i] <- covr_contr_rits[k-1, i] + ifelse(cond, 1, 0)
        width_contr_rits[k-1, i] <- width_contr_rits[k-1, i] + 
          abs(sim[[iter]]$contr[ate_ind[i], k-1, 3] - sim[[iter]]$contr[ate_ind[i], k-1, 2])
        bias_contr_rits[k-1, i] <- bias_contr_rits[k-1, i] + 
          sim[[iter]]$contr[ate_ind[i], k-1, 1] - contr_true[k-1]
        rmse_contr_rits[k-1, i] <- rmse_contr_rits[k-1, i] + 
          (sim[[iter]]$contr[ate_ind[i], k-1, 1] - contr_true[k-1])^2
      }
    }
    covr_ate_rits[k, i] <- covr_ate_rits[k, i] / n_iter
    width_ate_rits[k, i] <- width_ate_rits[k, i] / n_iter
    bias_ate_rits[k, i] <- bias_ate_rits[k, i] / n_iter
    
    if(k > 1){
      covr_contr_rits[k-1, i] <- covr_contr_rits[k-1, i] / n_iter
      width_contr_rits[k-1, i] <- width_contr_rits[k-1, i] / n_iter
      bias_contr_rits[k-1, i] <- bias_contr_rits[k-1, i] / n_iter
      rmse_contr_rits[k-1, i] <- sqrt(rmse_contr_rits[k-1, i] / n_iter)
    }
  }
}


# Bias + RMSE Table (requested by AE)
bias_tab <- matrix(NA, 4*(K-1), length(ind))
for(k in 1:(K-1)){
  bias_tab[4*(k-1)+1, ] <- bias_contr_std[k, ]
  bias_tab[4*(k-1)+2, ] <- bias_contr_rand[k, ]
  bias_tab[4*(k-1)+3, ] <- bias_contr_ts[k, ]
  bias_tab[4*(k-1)+4, ] <- bias_contr_rits[k, ]
}
rmse_tab <- matrix(NA, 4*(K-1), length(ind))
for(k in 1:(K-1)){
  rmse_tab[4*(k-1)+1, ] <- rmse_contr_std[k, ]
  rmse_tab[4*(k-1)+2, ] <- rmse_contr_rand[k, ]
  rmse_tab[4*(k-1)+3, ] <- rmse_contr_ts[k, ]
  rmse_tab[4*(k-1)+4, ] <- rmse_contr_rits[k, ]
}

est_err_tab <- matrix(paste(round(bias_tab, 2), "(", round(rmse_tab, 2), ")", sep = ""), 
       nrow = nrow(bias_tab))
colnames(est_err_tab) <- paste(ind)
rownames(est_err_tab) <- paste(c("Std", "Rand", "TS", "RiTS"), 
                               "(Arm", rep(2:4, each = 4), ")", sep = "")
xtable::xtable(est_err_tab)
# colnames(tab_metric) <- paste(ind)
# rownames(tab_metric) <- rep(c("Rand", "TS", "RiTS"), 2)
# print(xtable::xtable(tab_metric))

covr_contr <- cbind(covr_contr_rand, covr_contr_ts, covr_contr_rits)
rownames(covr_contr) <- paste("Arm", 2:K, "- Arm 1")
colnames(covr_contr) <- rep(paste(ind), 3)
print(xtable::xtable(covr_contr))


# width_contr <- cbind(width_contr_rand, width_contr_ts, width_contr_rits)
# rownames(width_contr) <- paste("Arm", 2:K, "- Arm 1")
# colnames(width_contr) <- rep(paste(ind), 3)
# print(xtable::xtable(width_contr))
# 
# 
# bias_contr <- cbind(bias_contr_rand, bias_contr_ts, bias_contr_rits)
# rownames(bias_contr) <- paste("Arm", 2:K, "- Arm 1")
# colnames(bias_contr) <- rep(paste(ind), 3)
# print(xtable::xtable(bias_contr))

# Proportion of times Arm 3 is selected for Phase III

n_comp <- 151
prop_top_std <- numeric(n_comp)
prop_top_rand <- numeric(n_comp)
prop_top_ts <- numeric(n_comp)
prop_top_rits <- numeric(n_comp)
for(i in 1:n_iter){
  prop_top_std <- prop_top_std + as.numeric(
    sapply(1:n_comp, function(t) which.max(rand_sim[[i]]$contr_standard[t, , 1]) == 3)
  )
  prop_top_rand <- prop_top_rand + as.numeric(
    sapply(1:n_comp, function(t) which.max(rand_sim[[i]]$contr[20+t, , 1]) == 3)
  )
  prop_top_ts <- prop_top_ts + as.numeric(
    sapply(1:n_comp, function(t) which.max(ts_sim[[i]]$contr[20+t, , 1]) == 3)
  )
  prop_top_rits <- prop_top_rits + as.numeric(
    sapply(1:n_comp, function(t) which.max(rits_sim[[i]]$contr[20+t, , 1]) == 3)
  )
}
prop_top_std <- prop_top_std / n_iter
prop_top_rand <- prop_top_rand / n_iter
prop_top_ts <- prop_top_ts / n_iter
prop_top_rits <- prop_top_rits / n_iter

plot(50:200, prop_top_std, type = "l")
lines(50:200, prop_top_rand, col = "red")
lines(50:200, prop_top_ts, col = "blue")
lines(50:200, prop_top_rits, col = "green")

zero_in_intv <- function(intv){
  a <- intv[1]; b <- intv[2]
  if(a < 0 && b > 0){
    return(TRUE)
  } else{
    return(FALSE)
  }
}

power_std <- numeric(n_comp)
power_rand <- numeric(n_comp)
power_ts <- numeric(n_comp)
power_rits <- numeric(n_comp)

for(i in 1:n_iter){
  power_std <- power_std + sapply(1:n_comp, function(t){
    rej <- sapply(1:3, function(k){
      zero_in_intv(rand_sim[[i]]$contr_standard[t, k, 2:3])
    })
    as.numeric(FALSE %in% rej)
  })
  power_rand <- power_rand + sapply(1:n_comp, function(t){
    rej <- sapply(1:3, function(k){
      zero_in_intv(rand_sim[[i]]$contr[t+20, k, 2:3])
    })
    as.numeric(FALSE %in% rej)
  })
  power_ts <- power_ts + sapply(1:n_comp, function(t){
    rej <- sapply(1:3, function(k){
      zero_in_intv(ts_sim[[i]]$contr[t+20, k, 2:3])
    })
    as.numeric(FALSE %in% rej)
  })
  power_rits <- power_rits + sapply(1:n_comp, function(t){
    rej <- sapply(1:3, function(k){
      zero_in_intv(rits_sim[[i]]$contr[t+20, k, 2:3])
    })
    as.numeric(FALSE %in% rej)
  })
}

power_std <- power_std / n_iter
power_rand <- power_rand / n_iter
power_ts <- power_ts / n_iter
power_rits <- power_rits / n_iter

plot(50:200, power_std, type = "l", ylim = c(0.95, 1))
lines(50:200, power_rand, col = "red")
lines(50:200, power_ts, col = "blue")
lines(50:200, power_rits, col = "green")


