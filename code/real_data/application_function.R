one_step_ts_batch_real <- function(x, t, param, eff, safe, K, d, log_dat,  
                                   min_prpn = 0.01, M = 1000, tr_start = 10, 
                                   tr_batch = 5, tr_lag = 10, ate_start = 20, 
                                   weight = 2, train = FALSE, rwd_sig = 0.1, 
                                   design = "no_clip"){
  # browser()
  last_tr_ind <- log_dat$last_tr_ind
  trt_hist <- log_dat$trt
  rwd_hist <- log_dat$reward
  cntx_hist <- log_dat$context
  prpns_mat <- log_dat$prpns_mat
  rwd_sig <- 1
  
  x <- matrix(x, ncol = 1)
  
  # online update
  beta_mean <- param$beta_mean
  beta_cov <- param$beta_cov
  beta_eff <- param$beta_eff
  if(t > tr_start){
    beta_mean[, , t] <- beta_mean[, , (t-1)]
    beta_cov[, , , t] <- beta_cov[, , , (t-1)]
  }
  param <- list(beta_mean = beta_mean, beta_cov = beta_cov, 
                beta_eff = beta_eff)
  # compute propensity scores
  prpns <- get_propensities_ts(t, x, K, d, param, tr_start = tr_start, M = M)
  
  if(design == "mad"){
    delt <- 1 / t^(0.24)
    prpns <- (1/K) * delt + (1 - delt) * prpns
  } else if(design == "clip"){
    prpns <- apply_floor(prpns, min_prpn)
  } else if(design == "mad_clip"){
    delt <- min_prpn * K
    prpns <- (1/K) * delt + (1 - delt) * prpns
  } else{ 
    prpns <- prpns # equivalent to "no_clip"
  }
  
  prpns_mat[t, ] <- prpns
  
  if(t <= 2*K){
    trt_new <- t %% K
    if(trt_new == 0) trt_new <- K
  } else{
    trt_new <- which.max(as.numeric(rmultinom(1, 1, prpns)))
  }
  
  
  # Generative Reward (Has to be modified)
  rwd_benf_new <- eff[t, trt_new]
  rwd_safe_new <- safe[t, trt_new]
  rwd_new <- rwd_benf_new
  
  # naive risk-adverse regret
  util <- (eff[t, ] + weight * safe[t, ]) / (weight+1)
  benf <- eff[t, ]
  saft <- safe[t, ]
  
  # Metrics
  regret_benf <- max(benf) - benf[trt_new]
  regret_safe <- max(saft) - saft[trt_new]
  regret <- max(util) - util[trt_new]
  
  subopt_cnt_benf <- as.numeric(which.max(benf) != trt_new)
  subopt_cnt_safe <- as.numeric(which.max(saft) != trt_new)
  subopt_cnt <- as.numeric(which.max(util) != trt_new)
  
  log_dat <- list(last_tr_ind = last_tr_ind,
                  trt = c(trt_hist, trt_new), 
                  reward = c(rwd_hist, rwd_new), 
                  context = rbind(cntx_hist, t(x)),
                  prpns_mat = prpns_mat)
  
  last_tr_ind <- log_dat$last_tr_ind
  trt_hist <- log_dat$trt
  rwd_hist <- log_dat$reward
  cntx_hist <- log_dat$context
  prpns_mat <- log_dat$prpns_mat
  
  if(train){
    # Updating Posterior Parameters
    for(i in (last_tr_ind+1):(t-tr_lag)){
      x <- as.numeric(cntx_hist[i, ])
      trt_batch <- trt_hist[i]
      rwd_batch <- rwd_hist[i]
      eff <- beta_eff[, , trt_batch]
      beta_eff[, , trt_batch] <- beta_eff[, , trt_batch] + 
        matrix(x, ncol = 1) %*% matrix(x, ncol = d) / rwd_sig^2
      beta_cov[, , trt_batch, t] <- solve(beta_eff[, , trt_batch])
      beta_mean[, trt_batch, t] <- beta_cov[, , trt_batch, t] %*% 
        (matrix(x, ncol = 1) * rwd_batch / rwd_sig^2 + 
           eff %*% beta_mean[, trt_batch, t])
    }
    last_tr_ind <- t-tr_lag
    log_dat["last_tr_dat"] <- last_tr_ind
  }
  
  # returns
  param <- list(beta_mean = beta_mean, beta_cov = beta_cov, 
                beta_eff = beta_eff)
  
  metrics <- list(regret = regret, regret_benf = regret_benf,
                  regret_safe = regret_safe,
                  subopt_cnt = subopt_cnt, subopt_cnt_benf = subopt_cnt_benf, 
                  subopt_cnt_safe = subopt_cnt_safe)
  
  list(trt_new = trt_new, rwd_new = rwd_new,
       param = param, metrics = metrics, log_dat = log_dat)
}


one_step_rits_batch_real <- function(x, t, param, eff, safe, K, d, 
                                log_dat, min_prpn = 0.01, M = 1000, 
                                tr_batch = 5, tr_start = 10, tr_lag = 10, 
                                ate_start = 20, weight = 2, train = FALSE, 
                                rwd_sig = 0.1, design = "no_clip"){
  last_tr_ind <- log_dat$last_tr_ind
  trt_hist <- log_dat$trt
  rwd_safe_hist <- log_dat$reward_safe
  rwd_benf_hist <- log_dat$reward_benf
  cntx_hist <- log_dat$context
  prpns_mat <- log_dat$prpns_mat
  rwd_sig <- 1
  
  x <- matrix(x, ncol = 1)
  # online update
  beta_mean <- param$beta_mean
  beta_cov <- param$beta_cov
  beta_eff <- param$beta_eff
  if(t > 1){
    beta_mean[, , t, ] <- beta_mean[, , (t-1), ]
    beta_cov[, , , t, ] <- beta_cov[, , , (t-1), ]
  }
  param <- list(beta_mean = beta_mean, beta_cov = beta_cov, 
                beta_eff = beta_eff)
  
  # compute propensity scores
  prpns <- get_propensities_rits(t, x, K, d, param, tr_start = tr_start, 
                                 weight = weight, last_tr_ind = last_tr_ind, 
                                 M = M)
  if(design == "mad"){
    delt <- 1 / t^(0.24)
    prpns <- (1/K) * delt + (1 - delt) * prpns
  } else if(design == "clip"){
    prpns <- apply_floor(prpns, min_prpn)
  } else if(design == "mad_clip"){
    delt <- min_prpn * K
    prpns <- (1/K) * delt + (1 - delt) * prpns
  } else{ 
    prpns <- prpns # equivalent to "no_clip"
  }
  prpns_mat[t, ] <- prpns
  
  if(t <= 2*K){
    trt_new <- t %% K
    if(trt_new == 0) trt_new <- K
  } else{
    trt_new <- which.max(as.numeric(rmultinom(1, 1, prpns)))
  }
  
  # Generative Reward (Has to be modified)
  rwd_benf_new <- eff[t, trt_new]
  rwd_safe_new <- safe[t, trt_new]
  
  # naive risk-adverse regret
  util <- (eff[t, ] + weight * safe[t, ]) / (weight+1)
  benf <- eff[t, ]
  saft <- safe[t, ]
  
  # Metrics
  regret_benf <- max(benf) - benf[trt_new]
  regret_safe <- max(saft) - saft[trt_new]
  regret <- max(util) - util[trt_new]
  
  subopt_cnt_benf <- as.numeric(which.max(benf) != trt_new)
  subopt_cnt_safe <- as.numeric(which.max(saft) != trt_new)
  subopt_cnt <- as.numeric(which.max(util) != trt_new)
  
  log_dat <- list(last_tr_ind = last_tr_ind,
                  trt = c(trt_hist, trt_new), 
                  reward_safe = c(rwd_safe_hist, rwd_safe_new),
                  reward_benf = c(rwd_benf_hist, rwd_benf_new),
                  context = rbind(cntx_hist, t(x)),
                  prpns_mat = prpns_mat)
  
  last_tr_ind <- log_dat$last_tr_ind
  trt_hist <- log_dat$trt
  rwd_safe_hist <- log_dat$reward_safe
  rwd_benf_hist <- log_dat$reward_benf
  cntx_hist <- log_dat$context
  prpns_mat <- log_dat$prpns_mat
  
  if(train){
    # browser()
    # Updating Posterior Parameters (Cont)
    for(e in 1:2){
      for(i in (last_tr_ind+1):(t-tr_lag)){
        x <- as.numeric(cntx_hist[i, ])
        trt_batch <- trt_hist[i]
        if(e == 1){
          rwd_batch <- rwd_benf_hist[i]
        } else{
          rwd_batch <- rwd_safe_hist[i]
        }
        eff <- beta_eff[, , trt_batch, e]
        beta_eff[, , trt_batch, e] <- beta_eff[, , trt_batch, e] + 
          matrix(x, ncol = 1) %*% matrix(x, ncol = d) / rwd_sig^2
        beta_cov[, , trt_batch, t, e] <- solve(beta_eff[, , trt_batch, e])
        beta_mean[, trt_batch, t, e] <- beta_cov[, , trt_batch, t, e] %*% 
          (matrix(x, ncol = 1) * rwd_batch / rwd_sig^2 + eff %*% 
             beta_mean[, trt_batch, t, e])
      }
    }
    last_tr_ind <- t
  }
  log_dat["last_tr_ind"] <- last_tr_ind
  
  # returns
  param <- list(beta_mean = beta_mean, beta_cov = beta_cov, 
                beta_eff = beta_eff)
  
  metrics <- list(regret = regret, regret_benf = regret_benf,
                  regret_safe = regret_safe,
                  subopt_cnt = subopt_cnt, subopt_cnt_benf = subopt_cnt_benf, 
                  subopt_cnt_safe = subopt_cnt_safe)
  
  list(trt_new = trt_new, rwd_safe_new = rwd_safe_new, 
       rwd_benf_new = rwd_benf_new, param = param, 
       metrics = metrics, log_dat = log_dat)
}



do_ts_batch_real <- function(X, eff, safe, K, tr_start = 30, tr_batch = 5,
                        weight = 1, tr_lag = 10, min_prpn = 0.01, ate_start = 30,
                        M = 1000, v = 10, placebo_arm = 1, alpha = 0.05, 
                        rwd_sig = 0.1, seed = 2024, design = "clip", 
                        first_peek = NULL){
  # browser()
  set.seed(seed)
  N <- nrow(X)
  d <- ncol(X)
  
  param <- list(beta_mean = array(0, c(d, K, N)), 
                beta_cov = array(0, c(d, d, K, N)), 
                beta_eff = array(0, c(d, d, K)))
  
  for(k in 1:K){
    param$beta_cov[, , k, 1] <- v * diag(d)
    param$beta_eff[, , k] <- (1/v) * diag(d)
  }
  log_dat <- list(last_tr_ind = 0, trt = numeric(0), 
                  reward = numeric(0), 
                  context = numeric(0),
                  prpns_mat = matrix(0, N, K))
  
  ate <- array(0, c(N, K, 3))
  contr <- array(0, c(N, K-1, 3))
  trt <- numeric(N)
  reward <- numeric(N)
  regret <- numeric(N)
  regret_benf <- numeric(N)
  regret_safe <- numeric(N)
  subopt_trt <- numeric(N)
  subopt_trt_benf <- numeric(N)
  subopt_trt_safe <- numeric(N)
  placebo_arm <- 1
  tr_first <- 0
  
  for(t in 1:N){
    if(t %% tr_batch == 0 && t >= tr_start){
      train <- TRUE
    } else{
      train <- FALSE
    }
    step <- one_step_ts_batch_real(X[t, ], t, param, eff, safe, K, d,
                              log_dat, tr_start = tr_start, tr_lag = tr_lag,
                              ate_start = ate_start, train = train, M = M, 
                              weight = weight, rwd_sig = rwd_sig, 
                              design = design, min_prpn = min_prpn)
    param <- step$param
    log_dat <- step$log_dat
    
    trt[t] <- step$trt_new
    reward[t] <- step$rwd_new
    regret[t] <- step$metrics$regret
    regret_benf[t] <- step$metrics$regret_benf
    regret_safe[t] <- step$metrics$regret_safe
    subopt_trt[t] <- step$metrics$subopt_cnt
    subopt_trt_benf[t] <- step$metrics$subopt_cnt_benf
    subopt_trt_safe[t] <- step$metrics$subopt_cnt_safe
  }
  times_ <- seq(ate_start, N)
  if(is.null(first_peek)) first_peek <- floor(N/2)
  asympcs <- get_asympcs(trt_hist = trt, rwd_hist = reward, 
                         prpns_mat = log_dat$prpns_mat, 
                         context_hist = X, placebo_arm = placebo_arm,
                         times = times_, alpha = alpha, first_peek = first_peek)
  ate <- asympcs[[1]]
  contr <- asympcs[[2]]
  
  list(trt = trt, reward = reward, regret = regret, regret_benf = regret_benf,
       regret_safe = regret_safe, subopt_trt = subopt_trt, 
       subopt_trt_benf = subopt_trt_benf, subopt_trt_safe = subopt_trt_safe,
       ate = ate, contr = contr, first_peek = first_peek, min_prpn = min_prpn,
       log_dat = log_dat, param = param, tr_first = tr_first)
}


do_rits_batch_real <- function(X, eff, safe, K, weight,
                          tr_start = 30, M = 1000, placebo_arm = 1,
                          tr_batch = 5, tr_lag = 10, ate_start = 30, 
                          v = 10, seed = 2024, design = "clip", alpha = 0.05,
                          rwd_sig = 0.1, min_prpn = 0.01, first_peek = NULL){
  # browser()
  set.seed(seed)
  N <- nrow(X)
  d <- ncol(X)
  E <- 2
  
  param <- list(beta_mean = array(0, c(d, K, N, E)), 
                beta_cov = array(0, c(d, d, K, N, E)), 
                beta_eff = array(0, c(d, d, K, E)))
  
  for(e in 1:E){
    for(k in 1:K){
      param$beta_cov[, , k, 1, e] <- v * diag(d)
      param$beta_eff[, , k, e] <- (1/v) * diag(d)
    }
  }
  log_dat <- list(last_tr_ind = 0, 
                  trt = numeric(0), 
                  reward_bin = numeric(0), 
                  reward_cont = numeric(0), 
                  context = numeric(0),
                  prpns_mat = matrix(0, N, K))
  
  ate <- array(0, c(N, K, 3))
  contr <- array(0, c(N, K-1, 3))
  trt <- numeric(N)
  reward_safe <- numeric(N)
  reward_benf <- numeric(N)
  regret <- numeric(N)
  regret_benf <- numeric(N)
  regret_safe <- numeric(N)
  subopt_trt <- numeric(N)
  subopt_trt_benf <- numeric(N)
  subopt_trt_safe <- numeric(N)
  placebo_arm <- 1
  tr_first <- 0
  
  for(t in 1:N){
    if(t %% tr_batch == 0 && t >= tr_start){
      train <- TRUE
    } else{
      train <- FALSE
    }
    step <- one_step_rits_batch_real(X[t, ], t, param, eff, safe, K, d, log_dat, 
                                tr_start = tr_start, tr_lag = tr_lag, 
                                tr_batch = tr_batch, M = M, 
                                ate_start = ate_start, train = train, 
                                weight = weight, rwd_sig = rwd_sig, 
                                design = design, min_prpn = min_prpn
                                )
    param <- step$param
    log_dat <- step$log_dat
    
    trt[t] <- step$trt_new
    reward_safe[t] <- step$rwd_safe_new
    reward_benf[t] <- step$rwd_benf_new
    regret[t] <- step$metrics$regret
    regret_benf[t] <- step$metrics$regret_benf
    regret_safe[t] <- step$metrics$regret_safe
    subopt_trt[t] <- step$metrics$subopt_cnt
    subopt_trt_benf[t] <- step$metrics$subopt_cnt_benf
    subopt_trt_safe[t] <- step$metrics$subopt_cnt_safe
  }
  
  times_ <- seq(ate_start, N)
  if(is.null(first_peek)) first_peek <- floor(N/2)
  asympcs <- get_asympcs(trt_hist = trt, rwd_hist = reward_benf, 
                         prpns_mat = log_dat$prpns_mat, 
                         context_hist = X, placebo_arm = placebo_arm,
                         times = times_, alpha = alpha, first_peek = first_peek)
  ate <- asympcs[[1]]
  contr <- asympcs[[2]]
  
  list(trt = trt, reward_safe = reward_safe, reward_benf = reward_benf, 
       regret = regret, regret_benf = regret_benf, regret_safe = regret_safe, 
       subopt_trt = subopt_trt, subopt_trt_benf = subopt_trt_benf, 
       subopt_trt_safe = subopt_trt_safe, ate = ate, contr = contr, 
       log_dat = log_dat, param = param, tr_first = tr_start)
}

get_metrics_rand <- function(eff, safe, arm, weight = 1){
  N <- length(arm)
  regret <- numeric(N); regret_benf <- regret; regret_safe <- regret
  subopt_cnt <- numeric(N); subopt_cnt_benf <- subopt_cnt; 
  subopt_cnt_safe <- subopt_cnt
  util <- (eff + safe * weight) / (weight + 1)
  for(i in 1:N){
    trt <- as.character(arm[i])
    regret[i] <- max(util[i, ]) - util[i, trt]
    regret_benf[i] <- max(eff[i, ]) - eff[i, trt]
    regret_safe[i] <- max(safe[i, ]) - safe[i, trt]
    trt_opt <- which.max(util[i, ])
    trt_opt_benf <- which.max(eff[i, ])
    trt_opt_safe <- which.max(safe[i, ])
    subopt_cnt[i] <- as.numeric(which.max(util[i, ]) != trt)
    subopt_cnt_benf[i] <- as.numeric(which.max(eff[i, ]) != trt)
    subopt_cnt_safe[i] <- as.numeric(which.max(safe[i, ]) != trt)
  }
  metrics <- list(trt = arm, regret = regret, regret_benf = regret_benf,
                  regret_safe = regret_safe,
                  subopt_trt = subopt_cnt, subopt_trt_benf = subopt_cnt_benf, 
                  subopt_trt_safe = subopt_cnt_safe)
  metrics
}

get_ate_rand <- function(obs_eff, arm, alpha = 0.05,
                         placebo_arm = 1, ate_start = 30){
  # browser()
  arm <- as.numeric(arm)
  N <- length(arm)
  K <- length(unique(arm))
  ate <- array(0, c(N, K, 3)); contr <- ate;
  avg <- matrix(0, N, K); std <- avg;
  for(k in 1:K){
    tmp <- sapply(ate_start:N, function(i){
      ind <- intersect(which(arm == k), 1:i)
      c(mean(obs_eff[ind]), sd(obs_eff[ind]))
    })
    avg[, k] <- c(rep(0, ate_start-1), tmp[1, ])
    std[, k] <- c(rep(0, ate_start-1), tmp[2, ])
    ate[, k, 1] <- avg[, k]
    width <- qnorm(1-alpha/2) * std[, k]
    ate[, k, 2] <- avg[, k] - width
    ate[, k, 3] <- avg[, k] + width
  }
  for(k in setdiff(1:K, placebo_arm)){
    contr[, k, 1] <- avg[, k] - avg[, placebo_arm]
    width <- qnorm(1-alpha/2) * sqrt(std[, k]^2 + std[, placebo_arm]^2)
    contr[, k, 2] <- avg[, k] - width
    contr[, k, 3] <- avg[, k] + width
  }
  contr <- contr[, -placebo_arm, ]
  list(ate = ate, contr = contr)
}
