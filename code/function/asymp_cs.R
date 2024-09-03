get_aipw_score <- function(t, K, trt_hist, rwd_hist, prpns_mat, tr_ind){
  if(t < 3 || length(tr_ind) == 0){
    stop("Too early to obtain AsympCS. Gather more data...")
  }
  ind1 <- tr_ind
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

get_aipw_var_running <- function(psi, tr_ind){
  t <- length(psi)
  var_running <- numeric(t)
  for(i in 1:t){
    var_running[i] <- sum((psi[1:i] - psi_running[i])^2) / i
  }
  var_running
}

get_rho <- function(first_peek, psi_var, alpha = 0.05){
  m <- first_peek; e <- exp(1)
  rho <- -2 * log(alpha) + log(-2 * log(alpha)) + 1
  rho <- rho / (psi_var * m * log(max(c(m, e))))
  rho <- sqrt(rho)
  rho
}

get_ugly_term <- function(t, rho, alpha){
  out <- 2 * (t * rho^2 + 1) * log(sqrt(t * rho^2 + 1) / alpha)
  out <- out / (t * rho)^2
  sqrt(out)
}

get_asympcs <- function(trt_hist, rwd_hist, prpns_mat,
                        trt_arm, placebo_arm, alpha = 0.05, first_peek = 50){
  # browser()
  T_ <- length(trt_hist); K <- length(unique(trt_hist)); m <- first_peek
  tr_ind <- sample(1:m, floor(m/2), replace = F)
  te_ind <- setdiff(1:m, tr_ind)
  e <- exp(1)
  asymp_cs <- matrix(0, T_, 3)
  
  for(t in m:T_){
    if(t > m){
      p <- rbinom(1, 1, 0.5)
      if(p == 0){
        tr_ind <- c(tr_ind, t)
      } else{
        te_ind <- c(te_ind, t)
      }
    }
    
    aipw_mat <- get_aipw_score(t, K, trt_hist, rwd_hist, prpns_mat, tr_ind)
    psi1 <- aipw_mat[tr_ind, trt_arm] - aipw_mat[tr_ind, placebo_arm]
    psi2 <- aipw_mat[te_ind, trt_arm] - aipw_mat[te_ind, placebo_arm]
    psi <- (sum(psi1) + sum(psi2))/t
    psi_var <- (var(psi1) + var(psi2))/2
    if(t == m){
      rho <- get_rho(first_peek, psi_var, alpha)
    }
    half_width <- sqrt(psi_var) * get_ugly_term(t, rho, alpha)
    asymp_cs[t, ] <- c(psi, psi - half_width, psi + half_width)
  }
  asymp_cs
}


# trt <- ts_out$trt; rwd <- ts_out$reward; prpns_mat <- ts_out$log_dat$prpns_mat
trt <- rits_out$trt; rwd <- rits_out$reward_benf; prpns_mat <- rits_out$log_dat$prpns_mat
par(mfrow = c(1, 3))
for(k in 2:4){
  tmp1 <- get_asympcs(trt, rwd, prpns_mat, k, 1, 0.05, first_peek = 30)
  
  plot(tmp1[, 1], type = "l", ylim = c(-3, 3), xlab = "Patient", ylab = "ATE")
  lines(tmp1[, 2])
  lines(tmp1[, 3])
  abline(h = mu_true[k] - mu_true[1])
}
par(mfrow = c(1, 1))

