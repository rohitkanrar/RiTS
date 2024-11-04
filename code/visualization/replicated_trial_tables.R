ts_sim <- readRDS("output/ts_sim_4.RData")
rand_sim <- readRDS("output/rand_sim_4.RData")
rits_sim <- readRDS("output/rits_sim_4.RData")
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
mu_true <- sim_dat$mu_true
# mu_true <- apply(X %*% beta_true[, , 1], 2, mean)
contr_true <- mu_true - mu_true[1]
contr_true <- contr_true[-1]

sim_choice <- readRDS("metadata/sim_choice.RData")
ate_start <- sim_choice$ate_start
ate_ind <- ind - ate_start

# RAND
sim <- rand_sim
covr_ate_rand <- matrix(0, K, length(ind))
width_ate_rand <- matrix(0, K, length(ind))
bias_ate_rand <- matrix(0, K, length(ind))
covr_contr_rand <- matrix(0, K-1, length(ind))
width_contr_rand <- matrix(0, K-1, length(ind))
bias_contr_rand <- matrix(0, K-1, length(ind))
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
        abs(sim[[iter]]$ate[ate_ind[i], k, 1] - mu_true[k])
      
      if(k > 1){
        cond <- (sim[[iter]]$contr[ate_ind[i], k-1, 2] <= contr_true[k-1]) &&
          (sim[[iter]]$contr[ate_ind[i], k-1, 3] >= contr_true[k-1])
        covr_contr_rand[k-1, i] <- covr_contr_rand[k-1, i] + ifelse(cond, 1, 0)
        width_contr_rand[k-1, i] <- width_contr_rand[k-1, i] + 
          abs(sim[[iter]]$contr[ate_ind[i], k-1, 3] - sim[[iter]]$contr[ate_ind[i], k-1, 2])
        bias_contr_rand[k-1, i] <- bias_contr_rand[k-1, i] + 
          abs(sim[[iter]]$contr[ate_ind[i], k-1, 1] - contr_true[k-1])
      }
    }
    covr_ate_rand[k, i] <- covr_ate_rand[k, i] / n_iter
    width_ate_rand[k, i] <- width_ate_rand[k, i] / n_iter
    bias_ate_rand[k, i] <- bias_ate_rand[k, i] / n_iter
    
    if(k > 1){
      covr_contr_rand[k-1, i] <- covr_contr_rand[k-1, i] / n_iter
      width_contr_rand[k-1, i] <- width_contr_rand[k-1, i] / n_iter
      bias_contr_rand[k-1, i] <- bias_contr_rand[k-1, i] / n_iter
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
        abs(sim[[iter]]$ate[ate_ind[i], k, 1] - mu_true[k])
      
      if(k > 1){
        cond <- (sim[[iter]]$contr[ate_ind[i], k-1, 2] <= contr_true[k-1]) &&
          (sim[[iter]]$contr[ate_ind[i], k-1, 3] >= contr_true[k-1])
        covr_contr_ts[k-1, i] <- covr_contr_ts[k-1, i] + ifelse(cond, 1, 0)
        width_contr_ts[k-1, i] <- width_contr_ts[k-1, i] + 
          abs(sim[[iter]]$contr[ate_ind[i], k-1, 3] - sim[[iter]]$contr[ate_ind[i], k-1, 2])
        bias_contr_ts[k-1, i] <- bias_contr_ts[k-1, i] + 
          abs(sim[[iter]]$contr[ate_ind[i], k-1, 1] - contr_true[k-1])
      }
    }
    covr_ate_ts[k, i] <- covr_ate_ts[k, i] / n_iter
    width_ate_ts[k, i] <- width_ate_ts[k, i] / n_iter
    bias_ate_ts[k, i] <- bias_ate_ts[k, i] / n_iter
    
    if(k > 1){
      covr_contr_ts[k-1, i] <- covr_contr_ts[k-1, i] / n_iter
      width_contr_ts[k-1, i] <- width_contr_ts[k-1, i] / n_iter
      bias_contr_ts[k-1, i] <- bias_contr_ts[k-1, i] / n_iter
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
        abs(sim[[iter]]$ate[ate_ind[i], k, 1] - mu_true[k])
      
      if(k > 1){
        cond <- (sim[[iter]]$contr[ate_ind[i], k-1, 2] <= contr_true[k-1]) &&
          (sim[[iter]]$contr[ate_ind[i], k-1, 3] >= contr_true[k-1])
        covr_contr_rits[k-1, i] <- covr_contr_rits[k-1, i] + ifelse(cond, 1, 0)
        width_contr_rits[k-1, i] <- width_contr_rits[k-1, i] + 
          abs(sim[[iter]]$contr[ate_ind[i], k-1, 3] - sim[[iter]]$contr[ate_ind[i], k-1, 2])
        bias_contr_rits[k-1, i] <- bias_contr_rits[k-1, i] + 
          abs(sim[[iter]]$contr[ate_ind[i], k-1, 1] - contr_true[k-1])
      }
    }
    covr_ate_rits[k, i] <- covr_ate_rits[k, i] / n_iter
    width_ate_rits[k, i] <- width_ate_rits[k, i] / n_iter
    bias_ate_rits[k, i] <- bias_ate_rits[k, i] / n_iter
    
    if(k > 1){
      covr_contr_rits[k-1, i] <- covr_contr_rits[k-1, i] / n_iter
      width_contr_rits[k-1, i] <- width_contr_rits[k-1, i] / n_iter
      bias_contr_rits[k-1, i] <- bias_contr_rits[k-1, i] / n_iter
    }
  }
}

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