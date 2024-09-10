get_propensities_ts <- function(t, x, K, d, param, M = 1000, tr_start = 10){
  if(t <= tr_start){
    return(rep(1/K, K))
  } else{
    prpns <- numeric(K)
    beta_mu <- param$beta_mean[, , t]
    beta_cov <- param$beta_cov[, , , t]
    samp <- matrix(0, M, K)
    for(k in 1:K){
      samp[, k] <- as.numeric(MASS::mvrnorm(M, beta_mu[, k], 
                                            beta_cov[, , k]) %*% x)
    }
    ind <- apply(samp, 1, which.max)
    for(k in 1:K){
      prpns[k] <- mean(ind == k)
    }
    return(prpns)
  }
}

get_propensities_rits <- function(t, x, K, d, param, tr_start = 10, 
                                  tr_lag = 10,  weight = 1, last_tr_ind = 0, 
                                  M = 1000){
  # if(t > tr_start) browser()
  if((t-tr_lag) < tr_start || last_tr_ind == 0){
    return(rep(1/K, K))
  } else{
    prpns <- numeric(K)
    samp_cont <- array(0, c(M, K, 2))
    samp <- matrix(0, M, K)
    beta_mu <- param$beta_mean[, , t, ]
    beta_cov <- param$beta_cov[, , , t, ]
    
    for(k in 1:K){
      for(e in 1:2){
        samp_cont[, k, e] <- as.numeric(MASS::mvrnorm(M, beta_mu[, k, e], 
                                                      beta_cov[, , k, e]) %*% x)
      }
      samp[, k] <- (samp_cont[, k, 2] * weight + 
                      samp_cont[, k, 1]) / (weight + 1)
    }
    
    ind <- apply(samp, 1, which.max)
    for(k in 1:K){
      prpns[k] <- mean(ind == k)
    }
    return(prpns)
  }
}

get_rwd_cont <- function(x, a, beta_true, sd = 0.1){
  avg_rwd <- as.numeric(crossprod(x, beta_true))
  r <- avg_rwd[a] + rnorm(1, sd = 0.1)
  list(avg_rwd = avg_rwd, rwd = r)
}

one_step_rand_biv <- function(x, K, beta_true, log_dat, weight, rwd_sig = 0.1){
  x <- matrix(as.numeric(x), ncol = 1)
  trt_hist <- log_dat$trt
  rwd_benf_hist <- log_dat$reward_benf
  rwd_safe_hist <- log_dat$reward_safe
  
  trt_new <- sample(1:K, 1)
  prpns <- rep(1/K, K)
  prpns_mat <- log_dat$prpns_mat
  prpns_mat <- rbind(prpns_mat, prpns)
  
  true_mu <- matrix(0, 2, K)
  tmp <- get_rwd_cont(x, trt_new, beta_true[, , 1])
  true_mu[1, ] <- tmp$avg_rwd
  rwd_benf_new <- tmp$rwd
  
  tmp <- get_rwd_cont(x, trt_new, beta_true[, , 2], rwd_sig)
  true_mu[2, ] <- tmp$avg_rwd
  rwd_safe_new <- tmp$rwd
  
  regret_benf <- max(true_mu[1, ]) - true_mu[1, trt_new]
  regret_safe <- max(true_mu[2, ]) - true_mu[2, trt_new]
  subopt_cnt_benf <- as.numeric(which.max(true_mu[1, ]) != trt_new)
  subopt_cnt_safe <- as.numeric(which.max(true_mu[2, ]) != trt_new)
  true_mu <- (true_mu[1, ] + weight * true_mu[2, ]) / (weight+1)
  
  regret <- max(true_mu) - true_mu[trt_new]
  subopt_cnt <- as.numeric(which.max(true_mu) != trt_new)
  
  metrics <- list(regret = regret, regret_safe = regret_safe,
                  regret_benf = regret_benf,
                  subopt_cnt = subopt_cnt, subopt_cnt_benf = subopt_cnt_benf,
                  subopt_cnt_safe = subopt_cnt_safe)
  log_dat <- list(trt = c(trt_hist, trt_new), 
                  reward_benf = c(rwd_benf_hist, rwd_benf_new),
                  reward_safe = c(rwd_safe_hist, rwd_safe_new),
                  prpns_mat = prpns_mat)
  
  list(trt_new = trt_new, rwd_benf_new = rwd_benf_new, 
       rwd_safe_new = rwd_safe_new,
       metrics = metrics, log_dat = log_dat)
}


one_step_ts_batch <- function(x, x_true, t, param, beta_true, K, d, log_dat,  
                              min_prpn = 0.01, M = 1000,
                              tr_start = 10, tr_lag = 10, ate_start = 20, 
                              weight = 2, train = FALSE, rwd_sig = 0.1,
                              design = "no_clip"){
  # browser()
  last_tr_ind <- log_dat$last_tr_ind
  trt_hist <- log_dat$trt
  rwd_hist <- log_dat$reward
  cntx_hist <- log_dat$context
  prpns_mat <- log_dat$prpns_mat
  
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
  # prpn_min <- floor_start / (t^floor_decay)
  # prpns <- apply_floor(prpns, prpn_min)
  if(design == "mad"){
    delt <- 1 / t^(0.24)
    prpns <- (1/K) * delt + (1 - delt) * prpns
  } else if(design == "clip"){
    prpns <- pmax(prpns, min_prpn)
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
  true_mu <- matrix(0, 2, K)
  tmp <- get_rwd_cont(x_true, trt_new, beta_true[, , 1], sd = rwd_sig)
  true_mu[1, ] <- tmp$avg_rwd
  rwd_benf_new <- tmp$rwd
  
  tmp <- get_rwd_cont(x_true, trt_new, beta_true[, , 2], sd = rwd_sig)
  true_mu[2, ] <- tmp$avg_rwd
  rwd_safe_new <- tmp$rwd
  rwd_new <- rwd_benf_new
  
  # naive risk-adverse regret
  true_mu_aug <- (true_mu[1, ] + weight * true_mu[2, ]) / (weight+1)
  
  # Metrics
  regret_benf <- max(true_mu[1, ]) - true_mu[1, trt_new]
  regret_safe <- max(true_mu[2, ]) - true_mu[2, trt_new]
  regret <- max(true_mu_aug) - true_mu_aug[trt_new]
  
  subopt_cnt_benf <- as.numeric(which.max(true_mu[1, ]) != trt_new)
  subopt_cnt_safe <- as.numeric(which.max(true_mu[2, ]) != trt_new)
  subopt_cnt <- as.numeric(which.max(true_mu_aug) != trt_new)
  
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


one_step_rits_batch <- function(x, x_true, t, param, beta_true, K, d, 
                                log_dat, floor_start = 0.005, 
                                min_prpn = 0.01, tr_start = 10, tr_lag = 10, 
                                M = 1000, ate_start = 20, tr_batch = 10, 
                                weight = 2, train = FALSE, rwd_sig = 0.1, 
                                design = "clip"){
  last_tr_ind <- log_dat$last_tr_ind
  trt_hist <- log_dat$trt
  rwd_safe_hist <- log_dat$reward_safe
  rwd_benf_hist <- log_dat$reward_benf
  cntx_hist <- log_dat$context
  prpns_mat <- log_dat$prpns_mat
  
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
                                 weight = weight, 
                                 last_tr_ind = last_tr_ind, M = M)
  # prpn_min <- floor_start / (t^floor_decay)
  # prpns <- apply_floor(prpns, prpn_min)
  if(design == "mad"){
    delt <- 1 / t^(0.24)
    prpns <- (1/K) * delt + (1 - delt) * prpns
  } else if(design == "clip"){
    prpns <- pmax(prpns, min_prpn)
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
  true_mu <- matrix(0, 2, K)
  for(e in 1:2){
    tmp <- get_rwd_cont(x_true, trt_new, beta_true[, , e], sd = rwd_sig)
    true_mu[e, ] <- tmp$avg_rwd
    if(e == 1) 
      rwd_benf_new <- tmp$rwd
    else
      rwd_safe_new <- tmp$rwd
  }
  
  
  # naive risk-adverse regret
  true_mu_aug <- (true_mu[1, ] + weight * true_mu[2, ]) / (weight+1)
  
  # Metrics
  regret_benf <- max(true_mu[1, ]) - true_mu[1, trt_new]
  regret_safe <- max(true_mu[2, ]) - true_mu[2, trt_new]
  regret <- max(true_mu_aug) - true_mu_aug[trt_new]
  
  subopt_cnt_benf <- as.numeric(which.max(true_mu[1, ]) != trt_new)
  subopt_cnt_safe <- as.numeric(which.max(true_mu[2, ]) != trt_new)
  subopt_cnt <- as.numeric(which.max(true_mu_aug) != trt_new)
  
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


do_rand_biv <- function(X_true, beta_true, tr_lag = 10, ate_start = 30, 
                        placebo_arm = 1, weight = 1, seed = 2024, 
                        rwd_sig = 0.1, alpha = 0.05, first_peek = 100){
  # browser()
  set.seed(seed)
  N <- nrow(X)
  K <- dim(beta_true)[2]
  log_dat <- list(trt = numeric(0), 
                  reward_benf = numeric(0),
                  reward_safe = numeric(0),
                  context = X, prpns_mat = numeric(0))
  
  trt <- numeric(N)
  reward_benf <- numeric(N)
  reward_safe <- numeric(N)
  regret <- numeric(N)
  regret_benf <- numeric(N)
  regret_safe <- numeric(N)
  subopt_trt <- numeric(N)
  subopt_trt_benf <- numeric(N)
  subopt_trt_safe <- numeric(N)
  
  for(t in 1:N){
    step <- one_step_rand_biv(X_true[t, ], K, beta_true, log_dat, weight, 
                              rwd_sig)
    log_dat <- step$log_dat
    trt[t] <- step$trt_new
    reward_benf[t] <- step$rwd_benf_new
    reward_safe[t] <- step$rwd_safe_new
    regret[t] <- step$metrics$regret
    regret_safe[t] <- step$metrics$regret_safe
    regret_benf[t] <- step$metrics$regret_benf
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
  
  list(trt = trt, reward_benf = reward_benf, reward_safe = reward_safe, 
       regret = regret, regret_safe = regret_safe, regret_benf = regret_benf, 
       subopt_trt = subopt_trt, subopt_trt_benf = subopt_trt_benf, 
       subopt_trt_safe = subopt_trt_safe, ate = ate, contr = contr,
       log_dat = log_dat)
}


do_ts_batch <- function(X, X_true, beta_true, tr_start = 20, tr_batch = 5, 
                        weight = 2, tr_lag = 10, ate_start = 30, M = 1000,
                        v = 10, placebo_arm = 1, min_prpn = 0.01, alpha = 0.05,
                        rwd_sig = 0.1, seed = 2024, design = "clip",
                        first_peek = NULL){
  # browser()
  set.seed(seed)
  N <- nrow(X)
  d <- ncol(X)
  K <- dim(beta_true)[2]
  
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
    step <- one_step_ts_batch(X[t, ], X_true[t, ], t, param, beta_true, K, d,
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
       regret_safe = regret_safe, subopt_trt = subopt_trt, alpha = alpha,
       subopt_trt_benf = subopt_trt_benf, subopt_trt_safe = subopt_trt_safe,
       ate = ate, contr = contr, first_peek = first_peek, min_prpn = min_prpn,
       log_dat = log_dat, param = param, tr_first = tr_first)
}

do_rits_batch <- function(X, X_true, beta_true, weight, rwd_sig = 0.1,
                          tr_start = 20, M = 1000, placebo_arm = 1,
                          tr_batch = 5, tr_lag = 10, ate_start = 30, 
                          v = 10, seed = 2024, design = "clip", alpha = 0.05,
                          min_prpn = 0.01, first_peek = NULL){
  # browser()
  set.seed(seed)
  N <- nrow(X)
  d <- ncol(X)
  K <- dim(beta_true)[2]
  E <- dim(beta_true)[3]
  
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
    step <- one_step_rits_batch(X[t, ], X_true[t, ], t, param, 
                                    beta_true, K, d, log_dat, 
                                    tr_start = tr_start, tr_lag = tr_lag, 
                                    tr_batch = tr_batch, M = M,
                                    ate_start = ate_start, 
                                    train = train, weight = weight, 
                                    rwd_sig = rwd_sig, design = design,
                                min_prpn = min_prpn)
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
       regret = regret, regret_benf = regret_benf,
       regret_safe = regret_safe, subopt_trt = subopt_trt, alpha = alpha,
       subopt_trt_benf = subopt_trt_benf, subopt_trt_safe = subopt_trt_safe, 
       ate = ate, contr = contr, first_peek = first_peek, min_prpn = min_prpn,
       log_dat = log_dat, param = param, tr_first = tr_first)
}