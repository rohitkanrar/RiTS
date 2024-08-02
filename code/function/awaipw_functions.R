# All the functions are translated to R from the python code available at
# https://github.com/gsbDBI/adaptive-confidence-intervals/blob/master/adaptive_CI/

apply_floor <- function(a, amin) {
  new <- pmax(a, amin)
  total_slack <- sum(new) - 1
  individual_slack <- new - amin
  c <- total_slack / sum(individual_slack)
  return(new - c * individual_slack)
}

twopoint_stable_var_ratio <- function(e, alpha){
  N <- nrow(e)
  K <- ncol(e)
  t <- matrix(1:N, ncol = 1)
  
  bad_lambda <- (1 - alpha) / ((1 - alpha) + N * (t / N) ^ alpha - t)
  good_lambda <- 1 / (1 + N - t)
  
  stopifnot(all(bad_lambda + 1e-7 >= good_lambda))
  
  lamb <- matrix(0, N, K)
  for(k in 1:K){
    lamb[, k] <- (1 - e[, k]) * bad_lambda + e[, k] * good_lambda
  }
  
  stopifnot(all(lamb >= 0))
  stopifnot(all(lamb <= 1 + 1e-8))
  for(k in 1:K){
    lamb[, k] <- pmax(0, pmin(lamb[, k], 1))
  }
  return(lamb)
}

stick_breaking <- function(Z){
  N <- nrow(Z)
  K <- ncol(Z)
  weights <- matrix(0, nrow = N, ncol = K)
  weight_sum <- rep(0, K)
  for (t in 1:N) {
    weights[t,] <- Z[t,] * (1 - weight_sum)
    weight_sum <- weight_sum + weights[t,]
  }
  return(weights)
}

evaluate_aipw_stats <- function(score, evalwts, alpha = 0.1) {
  estimate <- apply(evalwts * score, 2, sum) / apply(evalwts, 2, sum)
  stdr <- sqrt(apply(evalwts^2 * (score - estimate)^2, 2, sum)) / apply(evalwts, 2, sum)
  ci_radius <- qnorm(1 - alpha / 2) * stdr
  low <- estimate - ci_radius
  high <- estimate + ci_radius
  list(estimate = estimate, stderr = stdr, rad = ci_radius, 
       low = low, high = high)
}

get_arm_ate <- function(t, K, trt_hist, rwd_hist, prpns_mat, context, 
                        floor_decay = 0.7, alpha = 0.05){
  # browser()
  trt_hist <- trt_hist[1:t]
  rwd_hist <- rwd_hist[1:t]
  X <- context[1:t, ]
  prpns_mat <- prpns_mat[1:t, ]
  aipw_mat <- matrix(0, t, K)
  estm <- numeric(K)
  for(k in 1:K){
    ind <- which(trt_hist == k)
    estm[k] <- mean(rwd_hist[ind])
  }
  for(i in 1:t){
    for(k in 1:K){
      if(trt_hist[i] == k){
        aipw_mat[i, k] <- rwd_hist[i] / prpns_mat[i, k] +
          (1 - 1 / prpns_mat[i, k]) * estm[k]
      } else{
        aipw_mat[i, k] <- estm[k]
      }
    }
  }
  twopoint <- twopoint_stable_var_ratio(prpns_mat, floor_decay)
  twopoint_h2es <- stick_breaking(twopoint)
  wts_twopoint <- sqrt(twopoint_h2es * prpns_mat)
  wts_twopoint <- apply(wts_twopoint, 2, function(x) pmax(x, 0))
  
  evaluate_aipw_stats(aipw_mat, wts_twopoint, alpha)
}

# tmp <- get_arm_ate(100, K, log_dat_ts$trt, log_dat_ts$reward, 
#                    log_dat_ts$prpns_mat, log_dat_ts$context)

get_arm_contr <- function(arm_ate, placebo_arm = 1, alpha = 0.05){
  K <- length(arm_ate$estimate)
  contrs <- arm_ate$estimate - arm_ate$estimate[placebo_arm]
  contrs <- contrs[-placebo_arm]
  stdr <- arm_ate$stderr^2
  stdr <- stdr + stdr[placebo_arm]
  stdr <- sqrt(stdr[-placebo_arm])
  alpha <- 1 - alpha/(K-1)
  low <- contrs - qnorm(alpha) * stdr
  high <- contrs + qnorm(alpha) * stdr
  list(contrs = contrs, low = low, high = high)
}
# get_arm_ate(10, log_dat_ts$prpns_mat, log_dat_ts$aipw_mat)

get_arm_ate_rand <- function(trt, rwd, K, alpha = 0.05){
  ate <- numeric(K)
  stdr <- numeric(K)
  for(k in 1:K){
    ate[k] <- mean(rwd[trt == k])
    stdr[k] <- sd(rwd[trt == k])
  }
  list(estimate = ate, low = ate - qnorm(1-alpha/2) * stdr,
       high = ate + qnorm(1-alpha/2) * stdr, stderr = stdr)
}