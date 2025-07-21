require(drconfseq)
library(parallel)
basic_learner <- function(y, X, newX, ipw){
  # browser()
  # mod <- lm(y ~ ., data.frame(y = y, X = X))
  # mod <- glmnet::cv.glmnet(X, y, alpha = 0, nfolds = min(10, length(y)), grouped = F, 
  #                          nlambda = 20, weights = ipw)
  # tmp <- predict(mod, newx = newX, s = "lambda.1se")
  mod <- glmnet::glmnet(X, y, alpha = 0, nlambda = 10, lambda.max = 50, 
                        lambda.min.ratio = 0.2, weights = ipw)
  tmp <- predict(mod, newx = newX, s = 25)
  as.numeric(tmp)
}
main_effect_ridge <- function(y, X, newX, trt_ind, K, ipw){
  # browser()
  nte <- nrow(newX)
  trt_ind <- factor(trt_ind, levels = 1:K)
  X <- model.matrix(~0+trt_ind+X)
  X <- X[, -1]
  mod <- glmnet::glmnet(X, y, alpha = 0, weights = ipw)
  reg_est <- vector(mode = "list", length = K)
  tmp <- matrix(0, nte, K-1)
  newX <- cbind(tmp, newX)
  colnames(newX) <- colnames(X)
  for(k in 1:K){
    newX_k <- newX
    if(k != 1){
      newX_k[, k-1] <- 1
    }
    reg_est[[k]] <- as.numeric(predict(mod, newx = newX_k, s = 0.1))
  }
  reg_est
}

get_aipw_static <- function(y, reg_est, propensity, treatment, K){
  # browser()
  N <- length(y)
  aipw <- matrix(0, N, K)
  for(k in 1:K){
    trt_ind <- as.numeric(treatment == k)
    aipw[, k] <- reg_est[[k]] + (trt_ind / propensity[, k]) * (y - reg_est[[k]])
  }
  return(aipw)
}

# learner = c("full_ridge", "main_ridge", NULL)
# NULL correspond to IPW estimator.
get_aipw_seq <- function(treatment, y, propensity, X, 
                         train_idx = NULL, times = NULL, n_cores = 1,
                         cross_fit = TRUE, verbose = FALSE, learner = "full_ridge"){
  # browser()
  if (any(is.null(times))) {
    warning("\"times\" was left as null. Computing only at time n.")
    times = length(y)
  }
  K <- length(unique(treatment)); N <- length(y)
  if(!is.null(learner)){
    if (is.null(train_idx)){
      # train_idx <- rbinom(length(y), 1, 0.5)
      train_idx <- rep(c(0, 1), length.out = length(y))
    }
    train_idx <- train_idx == TRUE
    eval_idx <- 1 - train_idx == TRUE
    if (cross_fit) {
      train_indices <- list(train_idx, eval_idx)
    }
    else {
      train_indices <- list(train_idx)
    }
  }
  aipw_master <- vector(mode = "list", length = length(times))
  i <- 1
  for(time in times){
    # browser()
    aipw_master[[i]] <- matrix(NA, N, K)
    if(is.null(learner)){
      reg_est <- vector(mode = "list", length = K)
      for(k in 1:K){
        reg_est[[k]] <- numeric(length(y))
      }
      aipw_master[[i]] <- get_aipw_static(y = y, 
                                          reg_est = reg_est, 
                                          propensity = propensity, 
                                          treatment = treatment, K = K)
    } else{
      if (verbose) {
        print(paste("Fitting nuisance functions at time", 
                    time))
      }
      for(train_idx in train_indices){
        train_idx_t <- train_idx == TRUE
        train_idx_t[1:length(train_idx) > time] = FALSE
        eval_idx_t <- 1 - train_idx == TRUE
        eval_idx_t[1:length(train_idx) > time] = FALSE
        y_train <- y[train_idx_t]
        y_eval <- y[eval_idx_t]
        X_train <- X[train_idx_t, ]
        X_eval <- X[eval_idx_t, ]
        treatment_train <- treatment[train_idx_t]
        treatment_eval <- treatment[eval_idx_t]
        reg_est <- vector(mode = "list", length = K)
        prpn_train <- propensity[train_idx_t, ]
        prpn_eval <- propensity[eval_idx_t, ]
        if(learner == "main_ridge"){
          ipw <- length(prpn_train)
          for(k in 1:K){
            ipw[treatment_train == k] <- 1 / prpn_train[treatment_train == k, k]
          }
          reg_est <- main_effect_ridge(y = y_train, X = X_train, newX = X_eval, 
                                       trt_ind = treatment_train, K = K, ipw = ipw)
        } else if (learner == "full_ridge"){
          for(k in 1:K){
            # browser()
            reg_est[[k]] <- tryCatch({
              basic_learner(y = y_train[treatment_train == k], 
                            X = X_train[treatment_train == k, ], 
                            newX = X_eval, 
                            ipw = 1 / prpn_train[treatment_train == k, k])
            }, error = function(e){
              numeric(nrow(X_eval))
            })
          }
        } else{
          stop("learner should be one of c('main_ridge', 'full_ridge')")
        }
        # browser()
        aipw_master[[i]][eval_idx_t, ] <- get_aipw_static(y = y_eval, 
                                                          reg_est = reg_est, 
                                                          propensity = prpn_eval, 
                                                          treatment = treatment_eval, K = K)
      }
    }
    i <- i+1
  }
  names(aipw_master) <- times
  return(aipw_master)
}


get_asympcs <- function(trt_hist, rwd_hist, prpns_mat, context_hist,
                        placebo_arm, times, aipw_master = NULL,
                        alpha = 0.05, first_peek = 50, n_cores = 1, 
                        learner = "full_ridge"){
  # browser()
  if(times[1] > first_peek){
    times <- c(first_peek, times)
  }
  N <- length(trt_hist); K <- length(unique(trt_hist)); m <- first_peek
  asymp_cs <- matrix(0, N, 3)
  aipw_master <- get_aipw_seq(y = rwd_hist, X = context_hist, 
                              treatment = trt_hist, propensity = prpns_mat, 
                              times = times, learner = learner)
  rho2 <- drconfseq::best_rho2_exact(t_opt = m, alpha_opt = alpha)
  aipw_ate_list <- vector(mode = "list", length = length(aipw_master))
  aipw_contr_list <- aipw_ate_list
  ate_seq <- array(NA, c(length(times), K, 3))
  contr_seq <- array(NA, c(length(times), K-1, 3))
  k <- 1
  for(trt_arm in 1:K){
    i <- 1
    for(time in times){
      ate <- aipw_master[[paste(time)]][1:as.numeric(time), ]
      aipw_ate_list[[i]] <- ate[, trt_arm]
      ate <- ate[, trt_arm] - ate[, placebo_arm]
      aipw_contr_list[[i]] <- ate
      i <- i+1
    }
    names(aipw_ate_list) <- times
    names(aipw_contr_list) <- times
    
    asympcs_list <- mclapply(aipw_ate_list, function(aipw_ate) {
      acs <- drconfseq::lyapunov_asympcs(aipw_ate, rho2 = rho2, 
                                         alpha = alpha/K, # Bonferroni's correction
                                         return_all_times = FALSE)
      return(c((acs$l + acs$u)/2, acs$l, acs$u))
    }, mc.cores = n_cores)
    
    ate_seq[, k, ] <- do.call(rbind, asympcs_list)
    
    if(trt_arm == placebo_arm){
      k <- k+1
      next
    }
    asympcs_list <- mclapply(aipw_contr_list, function(aipw_contr) {
      acs <- drconfseq::lyapunov_asympcs(aipw_contr, rho2 = rho2, 
                                         alpha = alpha/(K-1), # Bonferroni's correction
                                         return_all_times = FALSE)
      return(c((acs$l + acs$u)/2, acs$l, acs$u))
    }, mc.cores = n_cores)
    
    contr_seq[, k-1, ] <- do.call(rbind, asympcs_list)
    k <- k+1
  }
  dimnames(ate_seq) <- list(
    Times = times, Arms = 1:K, CS = c("Center", "Lower", "Upper") 
  )
  dimnames(contr_seq) <- list(
    Times = times, Arms = sort(setdiff(1:K, placebo_arm)), 
    CS = c("Center", "Lower", "Upper") 
  )
  
  return(list(ate_seq, contr_seq))
}

add_asympcs <- function(out, ate_start, batch = 5, placebo_arm = 1, 
                        alpha = 0.05, first_peek = NULL, learner = "full_ridge"){
  if(is.null(first_peek)) first_peek <- ate_start
  trt <- out$trt; N <- length(trt)
  if(!is.null(out$reward_benf)){
    reward <- out$reward_benf
  } else{
    reward <- out$reward
  }
  times_ <- seq(ate_start, N, batch)
  asympcs <- get_asympcs(trt_hist = trt, rwd_hist = reward, 
                         prpns_mat = out$log_dat$prpns_mat, 
                         context_hist = out$log_dat$context, 
                         placebo_arm = placebo_arm, times = times_, 
                         alpha = alpha, first_peek = first_peek, 
                         learner = learner)
  ate <- asympcs[[1]]
  contr <- asympcs[[2]]
  out[["alpha"]] <- alpha
  if(is.null(learner)){
    out[["ate_ipw"]] <- ate
    out[["contr_ipw"]] <- contr
  } else{
    out[["ate"]] <- ate
    out[["contr"]] <- contr
  }
  if(is.null(out[["learner"]])){
    out[["learner"]] <- learner
  } else{
    out[["learner"]] <- c(out[["learner"]], learner)
  }
  out[["first_peek"]] <- first_peek
  out[["ate_start"]] <- ate_start
  return(out)
}

add_asympcs_sim <- function(out_list, ate_start, batch = 5, placebo_arm = 1, 
                            alpha = 0.05, first_peek = NULL, n_cores = 1,
                            force_compute = FALSE, learner = "full_ridge"){
  n_iter <- length(out_list)
  out_list <- mclapply(out_list, function(out){
    if(is.null(learner)){
      if(is.null(out$ate_ipw) || is.null(out$contr_ipw) || force_compute){
        out <- add_asympcs(out = out, ate_start = ate_start, batch = batch, 
                           placebo_arm = placebo_arm, alpha = alpha, 
                           first_peek = first_peek, learner = learner) 
      }
    } else{
      if(is.null(out$ate) || is.null(out$contr) || force_compute){
        out <- add_asympcs(out = out, ate_start = ate_start, batch = batch, 
                           placebo_arm = placebo_arm, alpha = alpha, 
                           first_peek = first_peek, learner = learner)
      }
    }
    return(out)
  }, mc.cores = n_cores)
}
# tmp_list1 <- add_asympcs_sim(tmp_list, 20)

# tmp <- get_asympcs(ts_out$trt, ts_out$reward, ts_out$log_dat$prpns_mat,
#                    ts_out$log_dat$context, 1, # times = seq(50, 1000, 10),
#                    times = 50:length(ts_out$trt),
#                    first_peek = 100)

# # # trt <- ts_out$trt; rwd <- ts_out$reward; prpns_mat <- ts_out$log_dat$prpns_mat
# true_mu <- get_true_avg_rwd(ts_out$log_dat$context, ts_out$trt, beta_true[, , 1])
# # # true_mu_running <- cumsum(true_mu) / (1:length(true_mu))
# # # trt <- ts_sim[[1]]$trt; rwd <- ts_sim[[1]]$reward; prpns_mat <- ts_sim[[1]]$log_dat$prpns_mat
# # trt <- rits_out$trt; rwd <- rits_out$reward_benf; prpns_mat <- rits_out$log_dat$prpns_mat
# par(mfrow = c(1, 3))
# for(k in 2:4){
#   tmp1 <- tmp[[k-1]]
# 
#   plot(tmp1[, 1], type = "l", ylim = c(-3, 3), xlab = "Patient", ylab = "ATE")
#   lines(tmp1[, 2])
#   # lines(tmp[, 3])
#   true_mu_running <- cumsum(true_mu[, k] - true_mu[, 1]) / (1:length(ts_out$trt))
#   lines(true_mu_running[51:length(ts_out$trt)], col = "red")
#   # abline(h = mu_true[k] - mu_true[1])
# }
# par(mfrow = c(1, 1))

