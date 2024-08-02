print(timestamp())
source("code/function/awaipw_functions.R")
source("code/function/main_function.R")
source("code/function/misc.R")

sim_choice <- readRDS("output/sim_choice.RData")
sim_dat <- readRDS("output/sim_dat.RData")

N <- sim_choice$N
K <- sim_choice$K # cannot be changed
d <- sim_choice$d # cannot be changed
E <- 2 # cannot be changed
X <- sim_dat$X_true
X_true <- sim_dat$X_true
beta_true <- sim_dat$beta_true
placebo_arm <- sim_dat$placebo_arm
weight <- sim_dat$weight
tr_start <- sim_choice$tr_start
ate_start <- sim_choice$ate_start
floor_start <- sim_choice$floor_start
floor_decay <- sim_choice$floor_decay
reward_sig <- sim_choice$reward_sig


n_iter <- 1000
ts_sim <- vector(mode = "list", length = n_iter)
rand_sim <- ts_sim
rits_sim <- ts_sim
ts_mis_sim <- ts_sim
rits_mis_sim <- ts_mis_sim

X <- X_true[, 1:2]

for(iter in 1:n_iter){
  print(iter)
  seed <- iter
  rand_sim[[iter]] <- do_rand_biv(X_true, beta_true, seed = seed, weight = weight, 
                                  placebo_arm = placebo_arm, rwd_sig = reward_sig)
  ts_sim[[iter]] <- do_ts_batch(X_true, X_true, beta_true, seed = seed, weight = weight, 
                                placebo_arm = placebo_arm, tr_start = tr_start, 
                                ate_start = ate_start, floor_start = floor_start, 
                                floor_decay = floor_decay, rwd_sig = reward_sig)
  rits_sim[[iter]] <- do_rits_batch(X_true, X_true, beta_true, weight = weight, seed = seed, 
                                        tr_start = tr_start, ate_start = ate_start,
                                        floor_start = floor_start, floor_decay = floor_decay,
                                        placebo_arm = placebo_arm, rwd_sig = reward_sig)
  
  ts_mis_sim[[iter]] <- do_ts_batch(X, X_true, beta_true, seed = seed, weight = weight, 
                                placebo_arm = placebo_arm, tr_start = tr_start, 
                                ate_start = ate_start, floor_start = floor_start, 
                                floor_decay = floor_decay, rwd_sig = reward_sig)
  rits_mis_sim[[iter]] <- do_rits_batch(X, X_true, beta_true, weight = weight, seed = seed, 
                                    tr_start = tr_start, ate_start = ate_start,
                                    floor_start = floor_start, floor_decay = floor_decay,
                                    placebo_arm = placebo_arm, rwd_sig = reward_sig)
  
}

saveRDS(ts_sim, "output/ts_sim.RData")
saveRDS(rand_sim, "output/rand_sim.RData")
saveRDS(rits_sim, "output/rits_sim.RData")
saveRDS(ts_mis_sim, "output/ts_sim.RData")
saveRDS(rits_mis_sim, "output/rits_sim.RData")
print(timestamp())
