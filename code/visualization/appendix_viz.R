source("code/visualization/functions.R")
### Section A.1 (Empirical Behavior of Cumulative Regret)
## Cumulative Regret Plot
# High SNR
ts_sim_high <- readRDS("output/ts_sim_high.RData")
rand_sim_high <- readRDS("output/rand_sim_high.RData")
rits_sim_high <- readRDS("output/rits_sim_high.RData")
n_iter <- length(ts_sim_high)
N <- length(ts_sim_high[[1]]$trt)
ind <- seq(20, 200, 20)
K <- length(unique(rand_sim_high[[1]]$trt))

sim_reg_plot1 <- gen_cum_reg_bwplot(out_rand = rand_sim_high,
                                    out_ts = ts_sim_high,
                                    out_rits = rits_sim_high,
                                    ind = ind, criteria = "Utility")
sim_reg_plot2 <- gen_cum_reg_bwplot(out_rand = rand_sim_high,
                                    out_ts = ts_sim_high,
                                    out_rits = rits_sim_high,
                                    ind = ind, criteria = "Efficacy")
sim_reg_plot3 <- gen_cum_reg_bwplot(out_rand = rand_sim_high,
                                    out_ts = ts_sim_high,
                                    out_rits = rits_sim_high,
                                    ind = ind, criteria = "Safety")
sim_regret_plot <- sim_reg_plot1 + sim_reg_plot2 + sim_reg_plot3 +
  plot_layout(ncol = 3, guides = "collect")
ggsave("plot/regret_sim_bwplot_high.jpg", height = 4, width = 12, units = "in")


# Low SNR
ts_sim_low <- readRDS("output/ts_sim_low.RData")
rand_sim_low <- readRDS("output/rand_sim_low.RData")
rits_sim_low <- readRDS("output/rits_sim_low.RData")

sim_reg_plot1 <- gen_cum_reg_bwplot(out_rand = rand_sim_low,
                                    out_ts = ts_sim_low,
                                    out_rits = rits_sim_low,
                                    ind = ind, criteria = "Utility")
sim_reg_plot2 <- gen_cum_reg_bwplot(out_rand = rand_sim_low,
                                    out_ts = ts_sim_low,
                                    out_rits = rits_sim_low,
                                    ind = ind, criteria = "Efficacy")
sim_reg_plot3 <- gen_cum_reg_bwplot(out_rand = rand_sim_low,
                                    out_ts = ts_sim_low,
                                    out_rits = rits_sim_low,
                                    ind = ind, criteria = "Safety")
sim_regret_plot <- sim_reg_plot1 + sim_reg_plot2 + sim_reg_plot3 +
  plot_layout(ncol = 3, guides = "collect")
ggsave("plot/regret_sim_bwplot_low.jpg", height = 4, width = 12, units = "in")


## Frequency of Arm Allocation Plot
# High SNR
alloc_plot1 <- gen_freq_arm_alloc(out_rand = rand_sim_high,
                                  out_ts = ts_sim_high,
                                  out_rits = rits_sim_high, 
                                  ylims = c(0, 150), titl = "High SNR")

# Low SNR
alloc_plot2 <- gen_freq_arm_alloc(out_rand = rand_sim_low,
                                  out_ts = ts_sim_low,
                                  out_rits = rits_sim_low, 
                                  ylims = c(0, 150), titl = "Low SNR")
alloc_plot <- alloc_plot1 + alloc_plot2 + 
  plot_layout(ncol = 2, guides = "collect")
ggsave("plot/freq_arm_alloc.jpg", height = 4, width = 10, units = "in")


### Section A.2 (Quality of AsympCS)
sim_choice <- readRDS("metadata/sim_choice.RData")
ate_start <- sim_choice$ate_start
ind <- c(seq(20, 50, 10), 70, 100, 150, 200)
ate_ind <- ind - ate_start + 1

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

cum_miscov_ts_high <- gen_cum_miscov_plot(out = ts_sim_high, ate_true = mu_true, 
                                          contr_true = contr_true, 
                                          alpha = alpha, ate_start = ate_start,
                                          titl = "TS")
cum_miscov_rits_high <- gen_cum_miscov_plot(out = rits_sim_high, ate_true = mu_true, 
                                          contr_true = contr_true, 
                                          alpha = alpha, ate_start = ate_start,
                                          titl = "RiTS")
cum_miscov_high <- cum_miscov_ts_high + cum_miscov_rits_high +
  plot_layout(ncol = 2, guides = "collect")
ggsave("plot/cum_miscov_contr_high.jpg", height = 4, width = 8, units = "in")


## Low SNR
mu_true <- sim_dat$mu_true
contr_true <- mu_true - mu_true[1]
contr_true <- contr_true[setdiff(1:K, sim_dat$placebo_arm)]

cum_miscov_ts_low <- gen_cum_miscov_plot(out = ts_sim_low, ate_true = mu_true, 
                                          contr_true = contr_true, 
                                          alpha = alpha, ate_start = ate_start,
                                          titl = "TS")
cum_miscov_rits_low <- gen_cum_miscov_plot(out = rits_sim_low, ate_true = mu_true, 
                                            contr_true = contr_true, 
                                            alpha = alpha, ate_start = ate_start,
                                            titl = "RiTS")
cum_miscov_low <- cum_miscov_ts_low + cum_miscov_rits_low +
  plot_layout(ncol = 2, guides = "collect")
ggsave("plot/cum_miscov_contr_low.jpg", height = 4, width = 8, units = "in")
