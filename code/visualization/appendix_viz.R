source("code/visualization/viz_function.R")
out_dir <- "output/"
# High SNR
ts_sim_high <- readRDS(paste(out_dir, "ts_sim_dgp_high_min_prpn_0.1_tr_start_24.RData", 
                             sep = ""))
rand_sim_high <- readRDS(paste(out_dir, "rand_sim_dgp_high_min_prpn_0.005_tr_start_24.RData", 
                               sep = ""))
rits_sim_high <- readRDS(paste(out_dir, "rits_sim_dgp_high_min_prpn_0.1_tr_start_24.RData", 
                               sep = ""))
# Low SNR
ts_sim_low <- readRDS(paste(out_dir, "ts_sim_dgp_low_min_prpn_0.1_tr_start_24.RData", 
                            sep = ""))
rand_sim_low <- readRDS(paste(out_dir, "rand_sim_dgp_low_min_prpn_0.005_tr_start_24.RData", 
                              sep = ""))
rits_sim_low <- readRDS(paste(out_dir, "rits_sim_dgp_low_min_prpn_0.1_tr_start_24.RData", 
                              sep = ""))
# Null Efficacy
ts_sim_null <- readRDS(paste(out_dir, "ts_sim_dgp_null_min_prpn_0.1_tr_start_24.RData", 
                             sep = ""))
rand_sim_null <- readRDS(paste(out_dir, "rand_sim_dgp_null_min_prpn_0.005_tr_start_24.RData", 
                               sep = ""))
rits_sim_null <- readRDS(paste(out_dir, "rits_sim_dgp_null_min_prpn_0.1_tr_start_24.RData", 
                               sep = ""))

sim_choice <- readRDS("metadata/sim_choice.RData")
n_iter <- length(ts_sim_high)
N <- length(ts_sim_high[[1]]$trt)
ind <- round(
  c(seq(50, sim_choice$N*0.625, 10), sim_choice$N*(3/4), sim_choice$N)
)
K <- length(unique(rand_sim_high[[1]]$trt))


### Section A.1 (Quality of AsympCS)
sim_choice <- readRDS("metadata/sim_choice.RData")
ate_start <- rand_sim_high[[1]]$ate_start
ate_ind <- sapply(ind, function(i){
  which(as.numeric(dimnames(rand_sim_high[[1]]$contr)[[1]]) == i)
})

## Box plot of Width for CS
df_high <- gen_width_df(out_rand = rand_sim_high, 
                        out_ts = ts_sim_high, out_rits = rits_sim_high, 
                        ate_ind = ate_ind)
df_low <- gen_width_df(out_rand = rand_sim_low, 
                        out_ts = ts_sim_low, out_rits = rits_sim_low, 
                        ate_ind = ate_ind)
df_null <- gen_width_df(out_rand = rand_sim_null, 
                       out_ts = ts_sim_null, out_rits = rits_sim_null, 
                       ate_ind = ate_ind)
sim_wid <- gen_width_bwplot(df_high = df_high, df_low = df_low, df_null = df_null,
                            ind = ind, ylims = list(c(0, 2.5), c(0, 1.5), c(0, 0.6)))
ggsave("plot/width_bwplot.jpg", plot = sim_wid, height = 6, width = 10, 
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

# Null Efficacy
mu_true <- sim_dat$mu_true
contr_true <- rep(0, length(mu_true) - 1)

df_null <- gen_bias_df(out_rand = rand_sim_null, 
                      out_ts = ts_sim_null, out_rits = rits_sim_null, 
                      ate_ind = ate_ind, contr_true = contr_true)

sim_bias <- gen_bias_bwplot(df_high = df_high, df_low = df_low, 
                            df_null = df_null, ind = ind,
                            ylims = list(c(-0.5, 0.5), c(-0.25, 0.25)))
ggsave("plot/bias_bwplot.jpg", plot = sim_bias, height = 6, width = 10, 
       units = "in")


### Section A.2 (Impact of Clipping)
source("code/function/misc.R")
alpha <- sim_choice$alpha
delay_aipw <- 57

## High SNR
mu_true <- sim_dat$mu_true * 2
contr_true <- mu_true - mu_true[1]
contr_true <- contr_true[setdiff(1:K, sim_dat$placebo_arm)]
df_high <- gen_cum_miscov_df(out_rand = rand_sim_high, out_ts = ts_sim_high, 
                             out_rits = rits_sim_high, mu_true = mu_true, 
                             contr_true = contr_true, delay_aipw = delay_aipw)
## Low SNR
mu_true <- sim_dat$mu_true
contr_true <- mu_true - mu_true[1]
contr_true <- contr_true[setdiff(1:K, sim_dat$placebo_arm)]
df_low <- gen_cum_miscov_df(out_rand = rand_sim_low, out_ts = ts_sim_low, 
                             out_rits = rits_sim_low, mu_true = mu_true, 
                             contr_true = contr_true, delay_aipw = delay_aipw)
## Null Efficacy
mu_true <- sim_dat$mu_true
contr_true <- rep(0, length(mu_true) - 1)
df_null <- gen_cum_miscov_df(out_rand = rand_sim_null, out_ts = ts_sim_null, 
                            out_rits = rits_sim_null, mu_true = mu_true, 
                            contr_true = contr_true, delay_aipw = delay_aipw)
cum_miscov <- gen_cum_miscov_plot(df_high = df_high, df_low = df_low, 
                                  df_null = df_null, alpha = sim_choice$alpha, 
                                  ate_start = sim_choice$ate_start)
ggsave("plot/cum_miscov_contr.jpg", cum_miscov, height = 9, width = 10, 
       units = "in")





rm(list = ls())
source("code/visualization/viz_function.R")
out_dir <- "output/"
### Section A.13 (Misspecified model)
# High SNR
ts_sim_high <- readRDS(paste(out_dir, "ts_mis_sim_dgp_high_min_prpn_0.1_tr_start_24.RData", 
                             sep = ""))
rand_sim_high <- readRDS(paste(out_dir, "rand_mis_sim_dgp_high_min_prpn_0.005_tr_start_24.RData", 
                               sep = ""))
rits_sim_high <- readRDS(paste(out_dir, "rits_mis_sim_dgp_high_min_prpn_0.1_tr_start_24.RData", 
                               sep = ""))
# Low SNR
ts_sim_low <- readRDS(paste(out_dir, "ts_mis_sim_dgp_low_min_prpn_0.1_tr_start_24.RData", 
                            sep = ""))
rand_sim_low <- readRDS(paste(out_dir, "rand_mis_sim_dgp_low_min_prpn_0.005_tr_start_24.RData", 
                              sep = ""))
rits_sim_low <- readRDS(paste(out_dir, "rits_mis_sim_dgp_low_min_prpn_0.1_tr_start_24.RData", 
                              sep = ""))

# Null Efficacy
ts_sim_null <- readRDS(paste(out_dir, "ts_mis_sim_dgp_null_min_prpn_0.1_tr_start_24.RData", 
                            sep = ""))
rand_sim_null <- readRDS(paste(out_dir, "rand_mis_sim_dgp_null_min_prpn_0.005_tr_start_24.RData", 
                              sep = ""))
rits_sim_null <- readRDS(paste(out_dir, "rits_mis_sim_dgp_null_min_prpn_0.1_tr_start_24.RData", 
                              sep = ""))

sim_choice <- readRDS("metadata/sim_choice.RData")
n_iter <- length(ts_sim_high)
N <- length(ts_sim_high[[1]]$trt)
ind <- round(
  c(seq(50, sim_choice$N*0.625, 10), sim_choice$N*(3/4), sim_choice$N)
)
K <- length(unique(rand_sim_high[[1]]$trt))

## Cumulative Regrets
criteria <- c("Utility", "Efficacy", "Safety")
df_high <- data.frame()
for(cr in criteria){
  df_high <- rbind(df_high, 
                   gen_cum_reg_df(out_rand = rand_sim_high, out_ts = ts_sim_high, 
                                  out_rits = rits_sim_high, ind = ind, criteria = cr)
  )
}
df_high["dgp"] <- "High-SNR"
df_low <- data.frame()
for(cr in criteria){
  df_low <- rbind(df_low, 
                  gen_cum_reg_df(out_rand = rand_sim_low, out_ts = ts_sim_low, 
                                 out_rits = rits_sim_low, ind = ind, criteria = cr)
  )
}
df_low["dgp"] <- "Low-SNR"
df_null <- data.frame()
for(cr in criteria){
  df_null <- rbind(df_null, 
                  gen_cum_reg_df(out_rand = rand_sim_null, out_ts = ts_sim_null, 
                                 out_rits = rits_sim_null, ind = ind, criteria = cr)
  )
}
df_null["dgp"] <- "Null-Efficacy"
sim_regret_plot <- gen_cum_reg_bwplot(df_high, df_low, df_null, ind)
ggsave("plot/regret_mis_sim_bwplot.jpg", height = 6, width = 10, units = "in")

## Frequency of Arm Allocation Plot
df_high <- gen_freq_arm_alloc_df(out_rand = rand_sim_high,
                                 out_ts = ts_sim_high,
                                 out_rits = rits_sim_high)
df_low <- gen_freq_arm_alloc_df(out_rand = rand_sim_low,
                                out_ts = ts_sim_low,
                                out_rits = rits_sim_low)
df_null <- gen_freq_arm_alloc_df(out_rand = rand_sim_null,
                                out_ts = ts_sim_null,
                                out_rits = rits_sim_null)
df_high[["dgp"]] <- "High-SNR"; df_low[["dgp"]] <- "Low-SNR"; df_null[["dgp"]] <- "Null-Efficacy"
alloc_plot <- gen_freq_arm_alloc(df_high = df_high, df_low = df_low, df_null = df_null,
                                 ylims = c(0, 150))
# Proportion of trials where Arm 4 is the winner (without confidence)
winner_high <- gen_winner_curve_df(rand_out = rand_sim_high, ts_out = ts_sim_high, 
                                   rits_out = rits_sim_high, true_best_arm = 4, include_ipw = FALSE)
winner_low <- gen_winner_curve_df(rand_out = rand_sim_low, ts_out = ts_sim_low, 
                                  rits_out = rits_sim_low, true_best_arm = 4, include_ipw = FALSE)
winner_null <- gen_winner_curve_df(rand_out = rand_sim_null, ts_out = ts_sim_null, 
                                  rits_out = rits_sim_null, true_best_arm = 4, include_ipw = FALSE)
winner_high[["dgp"]] <- "High-SNR"; winner_low[["dgp"]] <- "Low-SNR" 
winner_null[["dgp"]] <- "Null-Efficacy"
winner <- rbind(winner_high, winner_low, winner_null)
winner[["type"]] <- "Winner (Arm 4)"

# Proportion of trials where stopping criteria is met
power_high <- gen_power_curve_df(rand_out = rand_sim_high, ts_out = ts_sim_high, 
                                 rits_out = rits_sim_high, min_thresh = 0.1, include_ipw = FALSE)
power_low <- gen_power_curve_df(rand_out = rand_sim_low, ts_out = ts_sim_low, 
                                rits_out = rits_sim_low, min_thresh = 0.1, include_ipw = FALSE)
power_null <- gen_power_curve_df(rand_out = rand_sim_null, ts_out = ts_sim_null, 
                                rits_out = rits_sim_null, min_thresh = 0.1, include_ipw = FALSE)
power_high[["dgp"]] <- "High-SNR"; power_low[["dgp"]] <- "Low-SNR"
power_null[["dgp"]] <- "Null-Efficacy"
power_df <- rbind(power_high, power_low, power_null)
power_df[["type"]] <- "Stopping Criteria"

metric_plots <- gen_metrics_plot(df_winner = winner, df_power = power_df)
# ggsave("plot/metrics.jpg", height = 4, width = 6, units = "in")

metric_alloc_plot <- alloc_plot + metric_plots + plot_layout(ncol = 2)
ggsave("plot/metric_alloc_mis_plot.jpg", height = 4, width = 15, units = "in")


# Quality of AsympCS (misspecified)
sim_choice <- readRDS("metadata/sim_choice.RData")
ate_start <- rand_sim_high[[1]]$ate_start
ate_ind <- sapply(ind, function(i){
  which(as.numeric(dimnames(rand_sim_high[[1]]$contr)[[1]]) == i)
})

## Box plot of Width for CS
df_high <- gen_width_df(out_rand = rand_sim_high, 
                        out_ts = ts_sim_high, out_rits = rits_sim_high, 
                        ate_ind = ate_ind)
df_low <- gen_width_df(out_rand = rand_sim_low, 
                       out_ts = ts_sim_low, out_rits = rits_sim_low, 
                       ate_ind = ate_ind)
df_null <- gen_width_df(out_rand = rand_sim_null, 
                       out_ts = ts_sim_null, out_rits = rits_sim_null, 
                       ate_ind = ate_ind)
sim_wid <- gen_width_bwplot(df_high = df_high, df_low = df_low, df_null = df_null,
                            ind = ind, ylims = list(c(0, 3.25), c(0, 1.75), c(0, 0.6)))
ggsave("plot/width_mis_bwplot.jpg", plot = sim_wid, height = 6, width = 10, 
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

# Null Efficacy
mu_true <- sim_dat$mu_true
contr_true <- rep(0, length(mu_true) - 1)

df_null <- gen_bias_df(out_rand = rand_sim_null, 
                      out_ts = ts_sim_null, out_rits = rits_sim_null, 
                      ate_ind = ate_ind, contr_true = contr_true)

sim_bias <- gen_bias_bwplot(df_high = df_high, df_low = df_low, df_null = df_null,
                            ind = ind,
                            ylims = list(c(-0.75, 0.75), c(-0.38, 0.38)))
ggsave("plot/bias_mis_bwplot.jpg", plot = sim_bias, height = 6, width = 10, 
       units = "in")