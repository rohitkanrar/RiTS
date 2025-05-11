source("code/function/asymp_cs.R")
source("code/function/main_function.R")
source("code/function/misc.R")

set.seed(787)
N <- 200
K <- 4 # cannot be changed
d <- 3 # cannot be changed
E <- 2 # cannot be changed
x <- rnorm(N)
X_true <- cbind(1, x, x^2)
X <- X_true
beta_true <- array(0, c(d, K, E))
beta_true[, 1, 1] <- get_beta(1, 0.005, 0.5)
beta_true[, 2, 1] <- get_beta(1.35, 0.1, 0.5)
beta_true[, 3, 1] <- get_beta(1.35, 0.05, 0.5)
beta_true[, 4, 1] <- get_beta(1.6, 0.1, 0.5)

beta_true[, 1, 2] <- get_gamma(1, 0)
beta_true[, 2, 2] <- get_gamma(1, 0.005)
beta_true[, 3, 2] <- get_gamma(1, 0.05)
beta_true[, 4, 2] <- get_gamma(1, 0.3)
beta_true <- 5 * beta_true
(mu_true <- 5*c(1-1.25*0.005*2, 1.35-1.25*0.1*2, 1.35-1.25*0.05*2, 
              1.6-1.25*0.1*2))
(nu_true <- 5*c(1, 1-0.005, 1-0.05, 1-0.3))

beta_true <- (2/5) * beta_true
mu_true <- (2/5) * mu_true; nu_true <- (2/5) * nu_true

weight <- 1
(util_true <- (mu_true + nu_true * weight) / (1+weight))
placebo_arm <- 1
seed_ <- 1
tr_start <- 20
ate_start <- 20
min_prpn <- 0.05
first_peek <- 30
reward_sig <- 1
# correct
rand_out <- do_rand_biv(X = X_true, X_true = X_true, beta_true = beta_true, 
                                weight = weight, seed = seed_, 
                                rwd_sig = reward_sig, tr_start = tr_start,
                                asympcs = FALSE)
ts_out <- do_ts_batch(X = X_true, X_true = X_true, beta_true = beta_true, 
                              weight = weight, seed = seed_, rwd_sig = reward_sig,
                              tr_start = tr_start, tr_batch = 5, tr_lag = 10,
                              M = 1000, v = 10, min_prpn = min_prpn, asympcs = FALSE)
rits_out <- do_rits_batch(X = X_true, X_true = X_true, beta_true = beta_true, 
                                  weight = weight, seed = seed_, rwd_sig = reward_sig,
                                  tr_start = tr_start, tr_batch = 5, tr_lag = 10,
                                  M = 1000, v = 10, min_prpn = min_prpn, asympcs = FALSE)

# # miss
# X <- X_true[, 1:2]
# ts_mis_out <- do_ts_batch(X, X_true, beta_true, seed = seed_, weight = weight, 
#                       tr_start = tr_start, ate_start = ate_start,
#                       placebo_arm = placebo_arm,
#                       floor_start = floor_start, floor_decay = floor_decay,
#                       rwd_sig = reward_sig)
# rits_mis_out <- do_rits_batch(X, X_true, beta_true, weight = weight, seed = seed, 
#                           tr_start = tr_start, ate_start = ate_start,
#                           floor_start = floor_start, floor_decay = floor_decay,
#                           placebo_arm = placebo_arm, rwd_sig = reward_sig)
# rand_mis_out <- do_rand_biv(X_true, beta_true, seed = seed, weight = weight, 
#                         placebo_arm = placebo_arm)
# 
# sim_dat <- list(X = X, X_true = X_true, mu_true = mu_true, nu_true = nu_true,
#                 util_true = util_true, tr_start = tr_start,
#                 ate_start = ate_start, beta_true = beta_true,
#                 placebo_arm = placebo_arm, weight = weight)
# sim_choice <- list(N = N, K = K, d = d, tr_start = tr_start, 
#                    ate_start = ate_start, floor_start = floor_start,
#                    floor_decay = floor_decay, reward_sig = reward_sig)

# saveRDS(sim_choice, "output/sim_choice.RData")
# saveRDS(sim_dat, "output/sim_dat.RData")
# saveRDS(ts_out, "output/ts_out.RData")
# saveRDS(rits_out, "output/rits_out.RData")
# saveRDS(rand_out, "output/rand_out.RData")
# saveRDS(ts_mis_out, "output/ts_mis_out.RData")
# saveRDS(rits_mis_out, "output/rits_mis_out.RData")
# saveRDS(rand_mis_out, "output/rand_mis_out.RData")