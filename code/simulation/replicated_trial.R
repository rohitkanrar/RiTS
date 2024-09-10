print(timestamp())
source("code/function/asymp_cs.R")
source("code/function/main_function.R")
source("code/function/misc.R")

# out_dir <- "output/"
out_dir <- "/work/LAS/zhanruic-lab/rohitk/git_repos_data/RiTS/output/"

sim_choice <- readRDS("metadata/sim_choice.RData")
sim_dat <- readRDS("metadata/sim_dat.RData")

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
min_prpns <- sim_choice$min_prpns
design <- sim_choice$design
first_peek <- sim_choice$first_peek
reward_sig <- sim_choice$reward_sig
alpha <- sim_choice$alpha


n_iter <- 1000
ts_sim <- vector(mode = "list", length = n_iter)
rand_sim <- ts_sim
rits_sim <- ts_sim
# ts_mis_sim <- ts_sim
# rits_mis_sim <- ts_mis_sim

# X <- X_true[, 1:2]

i <- 1
for(min_prpn in min_prpns){
  for(iter in 1:n_iter){
    print(iter)
    seed <- iter
    if(i == 1){
      rand_sim[[iter]] <- do_rand_biv(X_true, beta_true, seed = seed, 
                                      weight = weight, 
                                      placebo_arm = placebo_arm, 
                                      rwd_sig = reward_sig, alpha = alpha,
                                      first_peek = first_peek)
    }
    
    ts_sim[[iter]] <- do_ts_batch(X_true, X_true, beta_true, seed = seed, weight = weight, 
                                  placebo_arm = placebo_arm, tr_start = tr_start, 
                                  ate_start = ate_start, rwd_sig = reward_sig, 
                                  design = design, min_prpn = min_prpns, 
                                  first_peek = first_peek, alpha = alpha)
    rits_sim[[iter]] <- do_rits_batch(X_true, X_true, beta_true, weight = weight, seed = seed, 
                                      tr_start = tr_start, ate_start = ate_start,
                                      placebo_arm = placebo_arm, rwd_sig = reward_sig, 
                                      design = design, min_prpn = min_prpns, 
                                      first_peek = first_peek, alpha = alpha)
    
    # ts_mis_sim[[iter]] <- do_ts_batch(X, X_true, beta_true, seed = seed, weight = weight, 
    #                               placebo_arm = placebo_arm, tr_start = tr_start, 
    #                               ate_start = ate_start, floor_start = floor_start, 
    #                               floor_decay = floor_decay, rwd_sig = reward_sig)
    # rits_mis_sim[[iter]] <- do_rits_batch(X, X_true, beta_true, weight = weight, seed = seed, 
    #                                   tr_start = tr_start, ate_start = ate_start,
    #                                   floor_start = floor_start, floor_decay = floor_decay,
    #                                   placebo_arm = placebo_arm, rwd_sig = reward_sig)
    
  }
  saveRDS(ts_sim, paste(out_dir, "ts_sim_", i, ".RData", sep = ""))
  saveRDS(rand_sim, paste(out_dir, "rand_sim_", i, ".RData", sep = ""))
  saveRDS(rits_sim, paste(out_dir, "rits_sim_", i, ".RData", sep = ""))
  # saveRDS(ts_mis_sim, "output/ts_mis_sim.RData")
  # saveRDS(rits_mis_sim, "output/rits_mis_sim.RData")
  i <- i + 1 
}

print(timestamp())
