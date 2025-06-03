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
ind <- c(ts_sim_high[[1]]$tr_first, seq(30, N/2, 15), 
         seq(N/2+25, N, 25))
K <- length(unique(rand_sim_high[[1]]$trt))


### Section A.2 (Quality of AsympCS)
sim_choice <- readRDS("metadata/sim_choice.RData")
ate_start <- rand_sim_high[[1]]$ate_start
ate_ind <- sapply(ind, function(i){
  which(as.numeric(dimnames(rand_sim_high[[1]]$contr)[[1]]) == i)
})

## Box plot of Width for CS
# High SNR
df_high <- gen_width_df(out_rand = rand_sim_high, 
                        out_ts = ts_sim_high, out_rits = rits_sim_high, 
                        ate_ind = ate_ind)
df_low <- gen_width_df(out_rand = rand_sim_low, 
                        out_ts = ts_sim_low, out_rits = rits_sim_low, 
                        ate_ind = ate_ind)
sim_wid <- gen_width_bwplot(df_high = df_high, df_low = df_low, 
                                 ind = ind)
ggsave("plot/width_bwplot.jpg", plot = sim_wid, height = 4, width = 12, 
       units = "in")


## Box plot for Bias
sim_dat <- readRDS("metadata/sim_dat.RData")
# High SNR
mu_true <- sim_dat$mu_true * 2
contr_true <- mu_true - mu_true[1]
contr_true <- contr_true[setdiff(1:K, sim_dat$placebo_arm)]

df_high <- gen_bias_df(out_rand = rand_sim_high, 
                        out_ts = ts_sim_high, out_rits = rits_sim_high, 
                        ate_ind = ate_ind, contr_true = contr_true)

# Low SNR
mu_true <- sim_dat$mu_true
contr_true <- mu_true - mu_true[1]
contr_true <- contr_true[setdiff(1:K, sim_dat$placebo_arm)]

df_low <- gen_bias_df(out_rand = rand_sim_low, 
                       out_ts = ts_sim_low, out_rits = rits_sim_low, 
                       ate_ind = ate_ind, contr_true = contr_true)

sim_bias <- gen_bias_bwplot(df_high = df_high, df_low = df_low, ind = ind)
ggsave("plot/bias_bwplot.jpg", plot = sim_bias, height = 4, width = 12, 
       units = "in")


### Section A.3 (Impact of Clipping)
source("code/function/misc.R")
alpha <- sim_choice$alpha

## High SNR
mu_true <- sim_dat$mu_true * 2
contr_true <- mu_true - mu_true[1]
contr_true <- contr_true[setdiff(1:K, sim_dat$placebo_arm)]
df_high <- gen_cum_miscov_df(out_rand = rand_sim_high, out_ts = ts_sim_high, 
                             out_rits = rits_sim_high, mu_true = mu_true, 
                             contr_true = contr_true)
## Low SNR
mu_true <- sim_dat$mu_true
contr_true <- mu_true - mu_true[1]
contr_true <- contr_true[setdiff(1:K, sim_dat$placebo_arm)]
df_low <- gen_cum_miscov_df(out_rand = rand_sim_low, out_ts = ts_sim_low, 
                             out_rits = rits_sim_low, mu_true = mu_true, 
                             contr_true = contr_true)
cum_miscov <- gen_cum_miscov_plot(df_high = df_high, df_low = df_low, 
                                  alpha = sim_choice$alpha, 
                                  ate_start = sim_choice$ate_start)
ggsave("plot/cum_miscov_contr.jpg", cum_miscov, height = 6, width = 12, 
       units = "in")


### Real Data Application
rand_real_sim <- readRDS("output/real_dat/rand_real_sim.RData")
ts_real_sim <- readRDS("output/real_dat/ts_real_sim.RData")
rits_real_sim <- readRDS("output/real_dat/rits_real_sim.RData")

n_iter <- length(rand_real_sim)
N <- length(rand_real_sim[[1]]$trt)
ind <- c(rand_real_sim[[1]]$ate_start, seq(rand_real_sim[[1]]$first_peek, N, 50))
K <- length(unique(rand_real_sim[[1]]$trt))


## Frequency of Arm Allocation Plot
df_real <- gen_freq_arm_alloc_df(out_rand = rand_real_sim,
                                 out_ts = ts_real_sim,
                                 out_rits = rits_real_sim)
alloc_plot <- gen_freq_arm_alloc(df_high = df_real, df_low = NULL, 
                                 ylims = c(0, 400))
# Proportion of trials where Arm 4 is the winner (without confidence)
winner <- gen_winner_curve_df(rand_out = rand_real_sim, ts_out = ts_real_sim, 
                                   rits_out = rits_real_sim, true_best_arm = 6, 
                              include_std = TRUE)
winner[["type"]] <- "Winner (Arm 6)"

# Proportion of trials where stopping criteria is met
power_df <- gen_power_curve_df(rand_out = rand_real_sim, ts_out = ts_real_sim, 
                                 rits_out = rits_real_sim, min_thresh = 0.1, 
                               include_std = TRUE)
power_df[["type"]] <- "Stopping Criteria"

metric_plots <- gen_metrics_plot(df_winner = winner, df_power = power_df, 
                                 dgp_exists = FALSE)
metric_alloc_plot <- alloc_plot + metric_plots + plot_layout(ncol = 2)
ggsave("plot/metric_alloc_plot_real.jpg", height = 4, width = 12, units = "in")