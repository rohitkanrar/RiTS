ts_sim <- readRDS("output/ts_sim.RData")
rand_sim <- readRDS("output/rand_sim.RData")
rits_sim <- readRDS("output/rits_sim.RData")
n_iter <- length(ts_sim)
N <- length(ts_sim[[1]]$trt)
ind <- c(seq(50, 100, 10), 150, 200)
K <- length(unique(rand_sim[[1]]$trt))

# average cumulative regret (both continuous end-points)
tab_metric <- matrix(0, 6, length(ind))

tmp <- t(sapply(1:n_iter, function(i){
  cumsum(rand_sim[[i]]$regret)[ind]
}))
tab_metric[1, ] <- apply(tmp, 2, mean)

tmp <- t(sapply(1:n_iter, function(i){
  cumsum(ts_sim[[i]]$regret)[ind]
}))
tab_metric[2, ] <- apply(tmp, 2, mean)

tmp <- t(sapply(1:n_iter, function(i){
  cumsum(rits_sim[[i]]$regret)[ind]
}))
tab_metric[3, ] <- apply(tmp, 2, mean)

tmp <- t(sapply(1:n_iter, function(i){
  cumsum(rand_sim[[i]]$subopt_trt)[ind]
}))
tab_metric[4, ] <- apply(tmp, 2, mean)

tmp <- t(sapply(1:n_iter, function(i){
  cumsum(ts_sim[[i]]$subopt_trt)[ind]
}))
tab_metric[5, ] <- apply(tmp, 2, mean)

tmp <- t(sapply(1:n_iter, function(i){
  cumsum(rits_sim[[i]]$subopt_trt)[ind]
}))
tab_metric[6, ] <- apply(tmp, 2, mean)

### Quality of Causal Effects

mu_true <- 5*c(1-1.25*0.005*2, 1.35-1.25*0.1*2, 1.35-1.25*0.05*2,
             1.6-1.25*0.1*2)
# mu_true <- apply(X %*% beta_true[, , 1], 2, mean)
contr_true <- mu_true - mu_true[1]
contr_true <- contr_true[-1]


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
    tmp <- which(sapply(1:n_iter, 
                        function(i) match(TRUE, sim[[i]]$ate[, k, 1] != 0)) <= ind[i])
    for(iter in tmp){
      cond <- (sim[[iter]]$ate[ind[i], k, 2] <= mu_true[k]) &&
        (sim[[iter]]$ate[ind[i], k, 3] >= mu_true[k])
      covr_ate_rand[k, i] <- covr_ate_rand[k, i] + ifelse(cond, 1, 0)
      width_ate_rand[k, i] <- width_ate_rand[k, i] + 
        abs(sim[[iter]]$ate[ind[i], k, 3] - sim[[iter]]$ate[ind[i], k, 2])
      bias_ate_rand[k, i] <- bias_ate_rand[k, i] + 
        abs(sim[[iter]]$ate[ind[i], k, 1] - mu_true[k])
      
      if(k > 1){
        cond <- (sim[[iter]]$contr[ind[i], k-1, 2] <= contr_true[k-1]) &&
          (sim[[iter]]$contr[ind[i], k-1, 3] >= contr_true[k-1])
        covr_contr_rand[k-1, i] <- covr_contr_rand[k-1, i] + ifelse(cond, 1, 0)
        width_contr_rand[k-1, i] <- width_contr_rand[k-1, i] + 
          abs(sim[[iter]]$contr[ind[i], k-1, 3] - sim[[iter]]$contr[ind[i], k-1, 2])
        bias_contr_rand[k-1, i] <- bias_contr_rand[k-1, i] + 
          abs(sim[[iter]]$contr[ind[i], k-1, 1] - contr_true[k-1])
      }
    }
    covr_ate_rand[k, i] <- covr_ate_rand[k, i] / length(tmp)
    width_ate_rand[k, i] <- width_ate_rand[k, i] / length(tmp)
    bias_ate_rand[k, i] <- bias_ate_rand[k, i] / length(tmp)
    
    if(k > 1){
      covr_contr_rand[k-1, i] <- covr_contr_rand[k-1, i] / length(tmp)
      width_contr_rand[k-1, i] <- width_contr_rand[k-1, i] / length(tmp)
      bias_contr_rand[k-1, i] <- bias_contr_rand[k-1, i] / length(tmp)
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
                        function(i) match(TRUE, sim[[i]]$ate[, k, 1] != 0)) <= ind[i])
    for(iter in tmp){
      cond <- (sim[[iter]]$ate[ind[i], k, 2] <= mu_true[k]) &&
        (sim[[iter]]$ate[ind[i], k, 3] >= mu_true[k])
      covr_ate_ts[k, i] <- covr_ate_ts[k, i] + ifelse(cond, 1, 0)
      width_ate_ts[k, i] <- width_ate_ts[k, i] + 
        abs(sim[[iter]]$ate[ind[i], k, 3] - sim[[iter]]$ate[ind[i], k, 2])
      bias_ate_ts[k, i] <- bias_ate_ts[k, i] + 
        abs(sim[[iter]]$ate[ind[i], k, 1] - mu_true[k])
      
      if(k > 1){
        cond <- (sim[[iter]]$contr[ind[i], k-1, 2] <= contr_true[k-1]) &&
          (sim[[iter]]$contr[ind[i], k-1, 3] >= contr_true[k-1])
        covr_contr_ts[k-1, i] <- covr_contr_ts[k-1, i] + ifelse(cond, 1, 0)
        width_contr_ts[k-1, i] <- width_contr_ts[k-1, i] + 
          abs(sim[[iter]]$contr[ind[i], k-1, 3] - sim[[iter]]$contr[ind[i], k-1, 2])
        bias_contr_ts[k-1, i] <- bias_contr_ts[k-1, i] + 
          abs(sim[[iter]]$contr[ind[i], k-1, 1] - contr_true[k-1])
      }
    }
    covr_ate_ts[k, i] <- covr_ate_ts[k, i] / length(tmp)
    width_ate_ts[k, i] <- width_ate_ts[k, i] / length(tmp)
    bias_ate_ts[k, i] <- bias_ate_ts[k, i] / length(tmp)
    
    if(k > 1){
      covr_contr_ts[k-1, i] <- covr_contr_ts[k-1, i] / length(tmp)
      width_contr_ts[k-1, i] <- width_contr_ts[k-1, i] / length(tmp)
      bias_contr_ts[k-1, i] <- bias_contr_ts[k-1, i] / length(tmp)
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
                        function(i) match(TRUE, sim[[i]]$ate[, k, 1] != 0)) <= ind[i])
    for(iter in tmp){
      cond <- (sim[[iter]]$ate[ind[i], k, 2] <= mu_true[k]) &&
        (sim[[iter]]$ate[ind[i], k, 3] >= mu_true[k])
      covr_ate_rits[k, i] <- covr_ate_rits[k, i] + ifelse(cond, 1, 0)
      width_ate_rits[k, i] <- width_ate_rits[k, i] + 
        abs(sim[[iter]]$ate[ind[i], k, 3] - sim[[iter]]$ate[ind[i], k, 2])
      bias_ate_rits[k, i] <- bias_ate_rits[k, i] + 
        abs(sim[[iter]]$ate[ind[i], k, 1] - mu_true[k])
      
      if(k > 1){
        cond <- (sim[[iter]]$contr[ind[i], k-1, 2] <= contr_true[k-1]) &&
          (sim[[iter]]$contr[ind[i], k-1, 3] >= contr_true[k-1])
        covr_contr_rits[k-1, i] <- covr_contr_rits[k-1, i] + ifelse(cond, 1, 0)
        width_contr_rits[k-1, i] <- width_contr_rits[k-1, i] + 
          abs(sim[[iter]]$contr[ind[i], k-1, 3] - sim[[iter]]$contr[ind[i], k-1, 2])
        bias_contr_rits[k-1, i] <- bias_contr_rits[k-1, i] + 
          abs(sim[[iter]]$contr[ind[i], k-1, 1] - contr_true[k-1])
      }
    }
    covr_ate_rits[k, i] <- covr_ate_rits[k, i] / length(tmp)
    width_ate_rits[k, i] <- width_ate_rits[k, i] / length(tmp)
    bias_ate_rits[k, i] <- bias_ate_rits[k, i] / length(tmp)
    
    if(k > 1){
      covr_contr_rits[k-1, i] <- covr_contr_rits[k-1, i] / length(tmp)
      width_contr_rits[k-1, i] <- width_contr_rits[k-1, i] / length(tmp)
      bias_contr_rits[k-1, i] <- bias_contr_rits[k-1, i] / length(tmp)
    }
  }
}