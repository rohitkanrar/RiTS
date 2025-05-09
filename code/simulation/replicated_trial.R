print(timestamp())
source("code/function/asymp_cs.R")
source("code/function/main_function.R")
source("code/function/misc.R")
source("code/function/classical_ci.R")

set.seed(673176)
out_dir <- "output/"

sim_choice <- readRDS("metadata/sim_choice.RData")
sim_dat <- readRDS("metadata/sim_dat.RData")
N <- sim_choice$N
K <- sim_choice$K # cannot be changed
d <- sim_choice$d # cannot be changed
E <- 2 # cannot be changed
X_true <- sim_dat$X_true
X <- X_true[, 1:2]
beta_true_low <- sim_dat$beta_true
beta_true_high <- 2 * beta_true_low
placebo_arm <- sim_dat$placebo_arm
weight <- sim_dat$weight
design <- sim_choice$design
reward_sig <- sim_choice$reward_sig

tr_starts <- sim_choice$tr_start
min_prpns <- sim_choice$min_prpns
dgps <- c("low", "high")
cases <- expand.grid(dgp = dgps, min_prpn = min_prpns, tr_start = tr_starts)
n_iter <- 1000

for(i in 1:nrow(cases)){
  dgp <- cases[i, "dgp"]; min_prpn <- cases[i, "min_prpn"]
  tr_start <- cases[i, "tr_start"]
  case_str <- paste("dgp", dgp, "min_prpn", min_prpn, "tr_start", tr_start, 
                    sep = "_")
  print(case_str)
  
  ts_sim <- vector(mode = "list", length = n_iter)
  rand_sim <- ts_sim; rits_sim <- ts_sim; ts_mis_sim <- ts_sim
  rits_mis_sim <- ts_mis_sim
  
  if(dgp == "low"){
    beta_true <- beta_true_low
  } else{
    beta_true <- beta_true_high
  }
  
  for(iter in 1:n_iter){
    # if(iter %% 100 == 0) print(iter)
    print(iter)
    rand_sim[[iter]] <- do_rand_biv(X_true = X_true, beta_true = beta_true, 
                                    weight = weight, seed = NULL, 
                                    rwd_sig = reward_sig, asympcs = FALSE)
    ts_sim[[iter]] <- do_ts_batch(X = X_true, X_true = X_true, beta_true = beta_true, 
                                  weight = weight, seed = NULL, rwd_sig = reward_sig,
                                  tr_start = tr_start, tr_batch = 5, tr_lag = 10,
                                  M = 1000, v = 10, min_prpn = min_prpn, asympcs = FALSE)
    rits_sim[[iter]] <- do_rits_batch(X = X_true, X_true = X_true, beta_true = beta_true, 
                                      weight = weight, seed = NULL, rwd_sig = reward_sig,
                                      tr_start = tr_start, tr_batch = 5, tr_lag = 10,
                                      M = 1000, v = 10, min_prpn = min_prpn, asympcs = FALSE)
    ts_mis_sim[[iter]] <- do_ts_batch(X = X, X_true = X_true, beta_true = beta_true, 
                                      weight = weight, seed = NULL, rwd_sig = reward_sig,
                                      tr_start = tr_start, tr_batch = 5, tr_lag = 10,
                                      M = 1000, v = 10, min_prpn = min_prpn, asympcs = FALSE)
    rits_mis_sim[[iter]] <- do_rits_batch(X = X, X_true = X_true, beta_true = beta_true, 
                                          weight = weight, seed = NULL, rwd_sig = reward_sig,
                                          tr_start = tr_start, tr_batch = 5, tr_lag = 10,
                                          M = 1000, v = 10, min_prpn = min_prpn, asympcs = FALSE)
  }
  # rand_sim <- add_standard_ci(out = rand_sim, times = seq(ate_start, N), 
  #                             placebo_arm = placebo_arm, alpha = alpha)
  saveRDS(ts_sim, paste(out_dir, "ts_sim_", case_str, ".RData", sep = ""))
  saveRDS(rand_sim, paste(out_dir, "rand_sim_", case_str, ".RData", sep = ""))
  saveRDS(rits_sim, paste(out_dir, "rits_sim_", case_str, ".RData", sep = ""))
  saveRDS(ts_mis_sim, paste(out_dir, "ts_mis_sim_", case_str, ".RData", sep = ""))
  saveRDS(rits_mis_sim, paste(out_dir, "rits_mis_sim_", case_str, ".RData", sep = ""))
}

print(timestamp())