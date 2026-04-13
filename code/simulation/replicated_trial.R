print(timestamp())
source("code/function/main_function.R")
source("code/function/misc.R")

sim_choice <- readRDS("metadata/sim_choice.RData")
sim_dat <- readRDS("metadata/sim_dat.RData")
dgps <- c("low", "high", "null")
tr_starts <- sim_choice$tr_start
min_prpns <- sim_choice$min_prpns
cases <- expand.grid(dgp = dgps, min_prpn = min_prpns, tr_start = tr_starts)
n_iter <- 1000

library(parallel)
num_cores <- 16
results <- mclapply(1:nrow(cases), function(i, cases, sim_choice, sim_dat, n_iter){
  # browser()
  out_dir <- "output/"
  N <- sim_choice$N
  K <- sim_choice$K # cannot be changed
  d <- sim_choice$d # cannot be changed
  E <- 2 # cannot be changed
  beta_true_low <- sim_dat$beta_true
  beta_true_high <- 2 * beta_true_low
  beta_true_null <- sim_dat$beta_true
  # Arms still have safety trend; All arms have placebo efficacy
  for(j in 1:K){
    if(j == 1) next
    beta_true_null[, j, 1] <- beta_true_null[, 1, 1]
  }
  placebo_arm <- sim_dat$placebo_arm
  weight <- sim_dat$weight
  design <- sim_choice$design
  reward_sig <- sim_choice$reward_sig
  
  dgp <- cases[i, "dgp"]; min_prpn <- cases[i, "min_prpn"]
  tr_start <- cases[i, "tr_start"]
  case_str <- paste("dgp", dgp, "min_prpn", min_prpn, "tr_start", tr_start, 
                    sep = "_")
  case_str0 <- paste("dgp", dgp, "min_prpn", 0.005, "tr_start", sim_choice$ate_start, 
                    sep = "_")
  rand_file_name <- paste(out_dir, "rand_sim_", case_str0, ".RData", sep = "")
  rand_mis_file_name <- paste(out_dir, "rand_mis_sim_", case_str0, ".RData", sep = "")
  
  ts_sim <- vector(mode = "list", length = n_iter)
  rand_sim <- ts_sim; rits_sim <- ts_sim; rand_mis_sim <- ts_sim; 
  ts_mis_sim <- ts_sim; rits_mis_sim <- ts_mis_sim
  if(dgp == "low"){
    beta_true <- beta_true_low
  } else if(dgp == "high"){
    beta_true <- beta_true_high
  } else{
    beta_true <- beta_true_null
  }
  for(iter in 1:n_iter){
    # if(iter %% 100 == 0) print(iter)
    seed_ <- iter
    set.seed(seed_)
    X_true <- get_X(N = N)
    X <- X_true[, 1:2]
    ts_sim[[iter]] <- do_ts_batch(X = X_true, X_true = X_true, beta_true = beta_true, 
                                  weight = weight, seed = seed_, rwd_sig = reward_sig,
                                  tr_start = tr_start, tr_batch = 5, tr_lag = 10,
                                  M = 1000, v = 10, min_prpn = min_prpn, asympcs = FALSE)
    rits_sim[[iter]] <- do_rits_batch(X = X_true, X_true = X_true, beta_true = beta_true, 
                                      weight = weight, seed = seed_, rwd_sig = reward_sig,
                                      tr_start = tr_start, tr_batch = 5, tr_lag = 10,
                                      M = 1000, v = 10, min_prpn = min_prpn, asympcs = FALSE)
    ts_mis_sim[[iter]] <- do_ts_batch(X = X, X_true = X_true, beta_true = beta_true, 
                                      weight = weight, seed = seed_, rwd_sig = reward_sig,
                                      tr_start = tr_start, tr_batch = 5, tr_lag = 10,
                                      M = 1000, v = 10, min_prpn = min_prpn, asympcs = FALSE)
    rits_mis_sim[[iter]] <- do_rits_batch(X = X, X_true = X_true, beta_true = beta_true, 
                                          weight = weight, seed = seed_, rwd_sig = reward_sig,
                                          tr_start = tr_start, tr_batch = 5, tr_lag = 10,
                                          M = 1000, v = 10, min_prpn = min_prpn, asympcs = FALSE)
    saveRDS(ts_sim, paste(out_dir, "ts_sim_", case_str, ".RData", sep = ""))
    saveRDS(rits_sim, paste(out_dir, "rits_sim_", case_str, ".RData", sep = ""))
    saveRDS(ts_mis_sim, paste(out_dir, "ts_mis_sim_", case_str, ".RData", sep = ""))
    saveRDS(rits_mis_sim, paste(out_dir, "rits_mis_sim_", case_str, ".RData", sep = ""))
    if(min_prpn == 0.005){
      rand_sim[[iter]] <- do_rand_biv(X = X_true, X_true = X_true, beta_true = beta_true, 
                                      weight = weight, seed = seed_, 
                                      rwd_sig = reward_sig, tr_start = tr_start,
                                      asympcs = FALSE)
      rand_mis_sim[[iter]] <- do_rand_biv(X = X, X_true = X_true, beta_true = beta_true, 
                                          weight = weight, seed = seed_, 
                                          rwd_sig = reward_sig, tr_start = tr_start,
                                          asympcs = FALSE)
    }
  }
  if(min_prpn == 0.005){
    saveRDS(rand_sim, paste(out_dir, "rand_sim_", case_str0, ".RData", sep = ""))
    saveRDS(rand_mis_sim, paste(out_dir, "rand_mis_sim_", case_str0, ".RData", sep = ""))
  }
  return(NULL)
}, mc.cores = num_cores, cases = cases, sim_dat = sim_dat, sim_choice = sim_choice, n_iter = n_iter)

print(timestamp())