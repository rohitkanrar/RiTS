source("code/visualization/functions.R")
### Section A.1 (Empirical Behavior of Cumulative Regret)
## Cumulative Regret Plot
# High SNR
ts_sim_high <- readRDS("output/ts_sim_dgp_high_min_prpn_0.05_tr_start_20.RData")
rand_sim_high <- readRDS("output/rand_sim_dgp_high_min_prpn_0.005_tr_start_20.RData")
rits_sim_high <- readRDS("output/rits_sim_dgp_high_min_prpn_0.05_tr_start_20.RData")
# Low SNR
ts_sim_low <- readRDS("output/ts_sim_dgp_low_min_prpn_0.05_tr_start_20.RData")
rand_sim_low <- readRDS("output/rand_sim_dgp_low_min_prpn_0.005_tr_start_20.RData")
rits_sim_low <- readRDS("output/rits_sim_dgp_low_min_prpn_0.05_tr_start_20.RData")

n_iter <- length(ts_sim_high)
N <- length(ts_sim_high[[1]]$trt)
ind <- seq(ts_sim_high[[1]]$tr_first, N, 20)
K <- length(unique(rand_sim_high[[1]]$trt))

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
sim_regret_plot <- gen_cum_reg_bwplot(df_high, df_low, ind)
ggsave("plot/regret_sim_bwplot.jpg", height = 4, width = 12, units = "in")


## Frequency of Arm Allocation Plot
df_high <- gen_freq_arm_alloc_df(out_rand = rand_sim_high,
                                 out_ts = ts_sim_high,
                                 out_rits = rits_sim_high)
df_low <- gen_freq_arm_alloc_df(out_rand = rand_sim_low,
                                out_ts = ts_sim_low,
                                out_rits = rits_sim_low)
df_high[["dgp"]] <- "High-SNR"; df_low[["dgp"]] <- "Low-SNR"
alloc_plot <- gen_freq_arm_alloc(df_high = df_high, df_low = df_low, 
                                 ylims = c(0, 70))
# Proportion of trials where Arm 4 is the winner (without confidence)
winner_high <- gen_winner_curve_df(rand_out = rand_sim_high, ts_out = ts_sim_high, 
                               rits_out = rits_sim_high, true_best_arm = 4)
winner_low <- gen_winner_curve_df(rand_out = rand_sim_low, ts_out = ts_sim_low, 
                               rits_out = rits_sim_low, true_best_arm = 4)
winner_high[["dgp"]] <- "High-SNR"; winner_low[["dgp"]] <- "Low-SNR"
winner <- rbind(winner_high, winner_low)
winner[["type"]] <- "Winner (Arm 4)"

# Proportion of trials where stopping criteria is met
power_high <- gen_power_curve_df(rand_out = rand_sim_high, ts_out = ts_sim_high, 
                                 rits_out = rits_sim_high, min_thresh = 0.1)
power_low <- gen_power_curve_df(rand_out = rand_sim_low, ts_out = ts_sim_low, 
                                 rits_out = rits_sim_low, min_thresh = 0.1)
power_high[["dgp"]] <- "High-SNR"; power_low[["dgp"]] <- "Low-SNR"
power_df <- rbind(power_high, power_low)
power_df[["type"]] <- "Stopping Criteria"

metric_plots <- gen_metrics_plot(df_winner = winner, df_power = power_df)
# ggsave("plot/metrics.jpg", height = 4, width = 6, units = "in")

metric_alloc_plot <- alloc_plot + metric_plots + plot_layout(ncol = 2)
ggsave("plot/metric_alloc_plot.jpg", height = 4, width = 12, units = "in")
