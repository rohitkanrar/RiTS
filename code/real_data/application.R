real_dat <- readRDS("data/clean_data.RData")
countfact_eff <- readRDS("output/real_dat/countfact_eff.RData")
countfact_safe <- readRDS("output/real_dat/countfact_safe.RData")
sim_choice <- readRDS("metadata/sim_choice.RData")

tr_start <- sim_choice$tr_start
ate_start <- sim_choice$ate_start
min_prpns <- sim_choice$min_prpns
design <- sim_choice$design
first_peek <- sim_choice$first_peek
reward_sig <- sim_choice$reward_sig
alpha <- sim_choice$alpha

X_real <- as.matrix(real_dat[c("baseline_salt", "baseline_lymph", "age", "bmi",
                               "albumin", "lymph_marker")])
X_real <- cbind(1, X_real, (as.numeric(real_dat$severity)-1))

source("code/function/main_function.R")
source("code/real_data/application_function.R")
source("code/function/misc.R")
source("code/function/asymp_cs.R")
set.seed(2024)

K <- length(unique(real_dat[["arm"]]))
v <- 10
weight <- 1
ts_real <- do_ts_batch_real(X = X_real, eff = countfact_eff, 
                            safe = countfact_safe, K = K, tr_start = tr_start,
                            tr_batch = 1, weight = weight, tr_lag = 0, 
                            ate_start = ate_start, M = 1000, placebo_arm = 1, 
                            v = v, min_prpn = min_prpns[3], alpha = alpha,
                            design = design, first_peek = first_peek, 
                            rwd_sig = reward_sig, seed = 100)
rits_real <- do_rits_batch_real(X = X_real, eff = countfact_eff, 
                            safe = countfact_safe, K = K, tr_start = tr_start,
                            tr_batch = 1, weight = weight, tr_lag = 0, 
                            ate_start = ate_start, M = 1000, placebo_arm = 1, 
                            v = v, min_prpn = min_prpns[3], alpha = alpha,
                            design = design, first_peek = first_peek, 
                            rwd_sig = reward_sig, seed = 100)
rand_real <- get_metrics_rand(countfact_eff, countfact_safe, real_dat[["arm"]])

saveRDS(ts_real, "output/ts_real.RData")
saveRDS(rits_real, "output/rits_real.RData")
saveRDS(rand_real, "output/rand_real.RData")

rm(list = ls())
