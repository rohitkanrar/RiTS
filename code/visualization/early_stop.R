require(ggplot2); require(reshape2); require(patchwork)
source("code/visualization/functions.R")
ts_sim_high <- readRDS("output/varyingX/ts_sim_dgp_high_min_prpn_0.05_tr_start_20.RData")
rand_sim_high <- readRDS("output/varyingX/rand_sim_dgp_high_min_prpn_0.005_tr_start_20.RData")
rits_sim_high <- readRDS("output/varyingX/rits_sim_dgp_high_min_prpn_0.05_tr_start_20.RData")

ts_sim_low <- readRDS("output/varyingX/ts_sim_dgp_low_min_prpn_0.05_tr_start_20.RData")
rand_sim_low <- readRDS("output/varyingX/rand_sim_dgp_low_min_prpn_0.005_tr_start_20.RData")
rits_sim_low <- readRDS("output/varyingX/rits_sim_dgp_low_min_prpn_0.05_tr_start_20.RData")

winner_high <- gen_winner_curve(rand_sim_high, ts_sim_high, rits_sim_high, 
                                titl = "High SNR", true_best_arm = 4)
winner_low <- gen_winner_curve(rand_sim_low, ts_sim_low, rits_sim_low, 
                               titl = "Low SNR", true_best_arm = 4)
winner_curve <- winner_high + winner_low + plot_layout(ncol = 2, guides = "collect")
ggsave("plot/prop_select.jpg", height = 4, width = 12, units = "in")


source("code/function/misc.R")
power_high <- gen_power_curve(rand_sim_high, ts_sim_high, rits_sim_high, 
                              titl = "High SNR", min_thresh = 0.1)
power_low <- gen_power_curve(rand_sim_low, ts_sim_low, rits_sim_low, 
                             titl = "Low SNR", min_thresh = 0.1)
power_curve <- power_high + power_low + plot_layout(ncol = 2, guides = "collect")
ggsave("plot/power.jpg", height = 4, width = 12, units = "in")

rm(list = ls())
