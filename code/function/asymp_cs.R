get_aipw_score <- function(t, K, trt_hist, rwd_hist, prpns_mat){
  ind1 <- sample(1:t, floor(t/2))
  ind2 <- setdiff(1:t, ind1)
  trt_hist1 <- trt_hist[ind1]
  trt_hist2 <- trt_hist[ind2]
  rwd_hist1 <- rwd_hist[ind1]
  rwd_hist2 <- rwd_hist[ind2]
  prpns_mat <- prpns_mat[1:t, ]
  aipw_mat <- matrix(NA, t, K)
  estm1 <- numeric(K); estm2 <- estm1
  for(k in 1:K){
    ind <- which(trt_hist1 == k)
    estm1[k] <- mean(rwd_hist1[ind])
    ind <- which(trt_hist2 == k)
    estm2[k] <- mean(rwd_hist2[ind])
  }
  for(i in 1:t){
    for(k in 1:K){
      estm <- ifelse(i %in% ind1, estm1[k], estm2[k])
      if(trt_hist[i] == k){
        aipw_mat[i, k] <- rwd_hist[i] / prpns_mat[i, k] +
          (1 - 1 / prpns_mat[i, k]) * estm
      } else{
        aipw_mat[i, k] <- estm
      }
    }
  }
  aipw_mat
}

get_aipw_var_running <- function(psi, psi_running){
  t <- length(psi)
  var_running <- numeric(t)
  for(i in 1:t){
    var_running[i] <- sum((psi[1:i] - psi_running[i])^2) / i
  }
  var_running
}

get_asympCS <- function(aipw_mat, trt_arm, placebo_arm, alpha,
                        first_peek = 50){
  t <- nrow(aipw_mat)
  m <- first_peek
  e <- exp(1)
  asymp_cs <- matrix(0, t, 3)
  psi <- aipw_mat[, trt_arm] - aipw_mat[, placebo_arm]
  psi_running <- cumsum(psi) / (1:t)
  var_running <- get_aipw_var_running(psi, psi_running)
  rho_m <- (-2 * log(alpha) + log(-2 * log(alpha)) + 1) / 
    (var_running[m] * m * log(max(m, e)))
  rho_m <- sqrt(rho_m)
  
  for(i in 1:t){
    if(i < first_peek) next
    # browser()
    tmp <- i * var_running[i] * rho_m^2 + 1
    tmp <- 2 * tmp * log(sqrt(tmp)/alpha) / (i * rho_m)^2
    tmp <- sqrt(tmp)
    asymp_cs[i, ] <- c(psi_running[i], psi_running[i] - tmp, 
                       psi_running[i] + tmp)
  }
  asymp_cs
}

# trt <- ts_out$trt; rwd <- ts_out$reward; prpns_mat <- ts_out$log_dat$prpns_mat
trt <- rits_out$trt; rwd <- rits_out$reward_benf; prpns_mat <- rits_out$log_dat$prpns_mat
par(mfrow = c(1, 3))
tmp <- get_aipw_score(length(trt), length(unique(trt)), trt, rwd, prpns_mat)
for(k in 2:4){
  tmp1 <- get_asympCS(tmp, k, 1, 0.05, first_peek = 30)
  
  plot(tmp1[, 1], type = "l", ylim = c(-3, 3), xlab = "Patient", ylab = "ATE")
  lines(tmp1[, 2])
  lines(tmp1[, 3])
  abline(h = mu_true[k] - mu_true[1])
}
par(mfrow = c(1, 1))

