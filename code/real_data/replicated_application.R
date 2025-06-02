print(timestamp())
source("code/function/main_function.R")
source("code/function/misc.R")
source("code/function/asymp_cs.R")
real_dat <- readRDS("data/clean_data.RData")
countfact_eff <- readRDS("output/real_dat/countfact_eff.RData")
countfact_safe <- readRDS("output/real_dat/countfact_safe.RData")
counterfact <- list(eff = countfact_eff, safe = countfact_safe)
X_real <- as.matrix(real_dat[c("baseline_salt", "baseline_lymph", "age", "bmi",
                               "albumin", "lymph_marker")])
X_real <- cbind(1, X_real, (as.numeric(real_dat$severity)-1))
n_iter <- 1000
ts_real_sim <- vector(mode = "list", length = n_iter)
rand_real_sim <- ts_real_sim; rits_real_sim <- ts_real_sim

library(parallel); num_cores <- 20
results <- mclapply(1:n_iter, function(i, X_real, counterfact){
  d <- ncol(X_real); N <- nrow(X_real); K <- ncol(countfact_eff)
  tr_start <- d*2*K; ate_start <- tr_start; min_prpn <- 0.05
  first_peek <- floor(ate_start * 1.5); reward_sig <- 1; alpha <- 0.05
  seed_ <- i
  set.seed(i)
  boot_ind <- sample(1:N, N, replace = TRUE)
  rand_real <- do_rand_biv(X = X_real[boot_ind, ], X_true = NULL, beta_true = NULL, 
                     weight = weight, seed = NULL, rwd_sig = reward_sig, 
                     tr_start = tr_start, asympcs = TRUE, ate_start = ate_start, 
                     first_peek = first_peek, setup = "real_data", 
                     counterfact = counterfact)
  ts_real <- do_ts_batch(X = X_real[boot_ind, ], X_true = NULL, beta_true = NULL, 
                         seed = NULL, rwd_sig = reward_sig,
                         tr_start = tr_start, tr_batch = 5, tr_lag = 10,
                         M = 1000, v = 10, min_prpn = 0.05, asympcs = TRUE, 
                         ate_start = ate_start, first_peek = first_peek, 
                         setup = "real_data", counterfact = counterfact)
  rits_real <- do_rits_batch(X = X_real[boot_ind, ], X_true = NULL, beta_true = NULL, 
                         weight = 1, seed = NULL, rwd_sig = reward_sig,
                         tr_start = tr_start, tr_batch = 5, tr_lag = 10,
                         M = 1000, v = 10, min_prpn = 0.05, asympcs = TRUE, 
                         ate_start = ate_start, first_peek = first_peek, 
                         setup = "real_data", counterfact = counterfact)
  list(rand_real = rand_real, ts_real = ts_real, rits_real = rits_real)
}, mc.cores = num_cores, X_real = X_real, counterfact = counterfact)

for(i in 1:n_iter){
  rand_real_sim[[i]] <- results[[i]]$rand_real
  ts_real_sim[[i]] <- results[[i]]$ts_real
  rits_real_sim[[i]] <- results[[i]]$rits_real
}

saveRDS(rand_real_sim, "output/real_dat/rand_real_sim.RData")
saveRDS(ts_real_sim, "output/real_dat/ts_real_sim.RData")
saveRDS(rits_real_sim, "output/real_dat/rits_real_sim.RData")
print(timestamp())