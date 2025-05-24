source("code/visualization/viz_function.R")
### Section A.1 (Empirical Behavior of Cumulative Regret)
## Cumulative Regret Plot
# High SNR
ts_sim_high <- readRDS("output/ts_sim_dgp_high_min_prpn_0.05_tr_start_24.RData")
rand_sim_high <- readRDS("output/rand_sim_dgp_high_min_prpn_0.005_tr_start_24.RData")
rits_sim_high <- readRDS("output/rits_sim_dgp_high_min_prpn_0.05_tr_start_24.RData")
# Low SNR
ts_sim_low <- readRDS("output/ts_sim_dgp_low_min_prpn_0.05_tr_start_24.RData")
rand_sim_low <- readRDS("output/rand_sim_dgp_low_min_prpn_0.005_tr_start_24.RData")
rits_sim_low <- readRDS("output/rits_sim_dgp_low_min_prpn_0.05_tr_start_24.RData")

n_iter <- length(ts_sim_high)
N <- length(ts_sim_high[[1]]$trt)
ind <- c(ts_sim_high[[1]]$tr_first, seq(30, N, 10))
K <- length(unique(rand_sim_high[[1]]$trt))


### Section A.2 (Quality of AsympCS)
sim_choice <- readRDS("metadata/sim_choice.RData")
ate_start <- rand_sim_high[[1]]$ate_start
ate_ind <- sapply(ind, function(i){
  which(as.numeric(dimnames(rand_sim_high[[1]]$contr)[[1]]) == i)
})

## Box plot of Width for CS
# High SNR
sim_wid_arm2 <- gen_width_bwplot(out_rand = rand_sim_high,
                                 out_ts = ts_sim_high,
                                 out_rits = rits_sim_high, 
                                 ate_ind = ate_ind, arm = 2, ylims = c(0, 15))
sim_wid_arm3 <- gen_width_bwplot(out_rand = rand_sim_high,
                                 out_ts = ts_sim_high,
                                 out_rits = rits_sim_high, 
                                 ate_ind = ate_ind, arm = 3, ylims = c(0, 15))
sim_wid_arm4 <- gen_width_bwplot(out_rand = rand_sim_high,
                                 out_ts = ts_sim_high,
                                 out_rits = rits_sim_high, 
                                 ate_ind = ate_ind, arm = 4, ylims = c(0, 15))
sim_wid_high <- sim_wid_arm2 + sim_wid_arm3 + sim_wid_arm4 +
  plot_layout(ncol = 3, guides = "collect")
ggsave("plot/width_bwplot_high.jpg", height = 4, width = 12, units = "in")

# Low SNR
sim_wid_arm2 <- gen_width_bwplot(out_rand = rand_sim_low,
                                 out_ts = ts_sim_low,
                                 out_rits = rits_sim_low, 
                                 ate_ind = ate_ind, arm = 2, ylims = c(0, 15))
sim_wid_arm3 <- gen_width_bwplot(out_rand = rand_sim_low,
                                 out_ts = ts_sim_low,
                                 out_rits = rits_sim_low, 
                                 ate_ind = ate_ind, arm = 3, ylims = c(0, 15))
sim_wid_arm4 <- gen_width_bwplot(out_rand = rand_sim_low,
                                 out_ts = ts_sim_low,
                                 out_rits = rits_sim_low, 
                                 ate_ind = ate_ind, arm = 4, ylims = c(0, 15))
sim_wid_low <- sim_wid_arm2 + sim_wid_arm3 + sim_wid_arm4 +
  plot_layout(ncol = 3, guides = "collect")
ggsave("plot/width_bwplot_low.jpg", height = 4, width = 12, units = "in")


## Box plot for Bias
sim_dat <- readRDS("metadata/sim_dat.RData")

# High SNR
mu_true <- sim_dat$mu_true * 2
contr_true <- mu_true - mu_true[1]
contr_true <- contr_true[setdiff(1:K, sim_dat$placebo_arm)]

sim_bias_arm2 <- gen_bias_bwplot(out_rand = rand_sim_high,
                                 out_ts = ts_sim_high,
                                 out_rits = rits_sim_high, 
                                 ate_ind = ate_ind, arm = 2, 
                                 contr_true = contr_true, ylims = c(-2, 2))
sim_bias_arm3 <- gen_bias_bwplot(out_rand = rand_sim_high,
                                 out_ts = ts_sim_high,
                                 out_rits = rits_sim_high, 
                                 ate_ind = ate_ind, arm = 3, 
                                 contr_true = contr_true, ylims = c(-2, 2))
sim_bias_arm4 <- gen_bias_bwplot(out_rand = rand_sim_high,
                                 out_ts = ts_sim_high,
                                 out_rits = rits_sim_high, 
                                 ate_ind = ate_ind, arm = 4, 
                                 contr_true = contr_true, ylims = c(-2, 2))
sim_bias_high <- sim_bias_arm2 + sim_bias_arm3 + sim_bias_arm4 +
  plot_layout(ncol = 3, guides = "collect")
ggsave("plot/bias_bwplot_high.jpg", height = 4, width = 12, units = "in")

# Low SNR
mu_true <- sim_dat$mu_true
contr_true <- mu_true - mu_true[1]
contr_true <- contr_true[setdiff(1:K, sim_dat$placebo_arm)]

sim_bias_arm2 <- gen_bias_bwplot(out_rand = rand_sim_low,
                                 out_ts = ts_sim_low,
                                 out_rits = rits_sim_low, 
                                 ate_ind = ate_ind, arm = 2, 
                                 contr_true = contr_true, ylims = c(-2, 2))
sim_bias_arm3 <- gen_bias_bwplot(out_rand = rand_sim_low,
                                 out_ts = ts_sim_low,
                                 out_rits = rits_sim_low, 
                                 ate_ind = ate_ind, arm = 3, 
                                 contr_true = contr_true, ylims = c(-2, 2))
sim_bias_arm4 <- gen_bias_bwplot(out_rand = rand_sim_low,
                                 out_ts = ts_sim_low,
                                 out_rits = rits_sim_low, 
                                 ate_ind = ate_ind, arm = 4, 
                                 contr_true = contr_true, ylims = c(-2, 2))
sim_bias_low <- sim_bias_arm2 + sim_bias_arm3 + sim_bias_arm4 +
  plot_layout(ncol = 3, guides = "collect")
ggsave("plot/bias_bwplot_low.jpg", height = 4, width = 12, units = "in")


### Section A.3 (Impact of Clipping)
source("code/function/misc.R")
library(reshape2)
# To be done later after running the simulation for all values of clipping param.
# For now only one clipping parameter is considered.
alpha <- sim_choice$alpha

## High SNR
mu_true <- sim_dat$mu_true * 2
contr_true <- mu_true - mu_true[1]
contr_true <- contr_true[setdiff(1:K, sim_dat$placebo_arm)]

cum_miscov_rand_high <- gen_cum_miscov_plot(out = rand_sim_high, ate_true = mu_true, 
                                            contr_true = contr_true, 
                                            alpha = alpha, ate_start = ate_start,
                                            titl = "Rand")
cum_miscov_ts_high <- gen_cum_miscov_plot(out = ts_sim_high, ate_true = mu_true, 
                                          contr_true = contr_true, 
                                          alpha = alpha, ate_start = ate_start,
                                          titl = "TS")
cum_miscov_rits_high <- gen_cum_miscov_plot(out = rits_sim_high, ate_true = mu_true, 
                                            contr_true = contr_true, 
                                            alpha = alpha, ate_start = ate_start,
                                            titl = "RiTS")
cum_miscov_high <- cum_miscov_rand_high + cum_miscov_ts_high + 
  cum_miscov_rits_high + plot_layout(ncol = 3, guides = "collect")
ggsave("plot/cum_miscov_contr_high.jpg", height = 4, width = 12, units = "in")


## Low SNR
mu_true <- sim_dat$mu_true
contr_true <- mu_true - mu_true[1]
contr_true <- contr_true[setdiff(1:K, sim_dat$placebo_arm)]

cum_miscov_rand_low <- gen_cum_miscov_plot(out = rand_sim_low, ate_true = mu_true, 
                                           contr_true = contr_true, 
                                           alpha = alpha, ate_start = ate_start,
                                           titl = "Rand")
cum_miscov_ts_low <- gen_cum_miscov_plot(out = ts_sim_low, ate_true = mu_true, 
                                         contr_true = contr_true, 
                                         alpha = alpha, ate_start = ate_start,
                                         titl = "TS")
cum_miscov_rits_low <- gen_cum_miscov_plot(out = rits_sim_low, ate_true = mu_true, 
                                           contr_true = contr_true, 
                                           alpha = alpha, ate_start = ate_start,
                                           titl = "RiTS")
cum_miscov_low <- cum_miscov_rand_low + cum_miscov_ts_low + cum_miscov_rits_low +
  plot_layout(ncol = 3, guides = "collect")
ggsave("plot/cum_miscov_contr_low.jpg", height = 4, width = 12, units = "in")
