source("code/visualization/viz_function.R")
out_dir <- "output/"
# High SNR
ts_sim_high <- readRDS(paste(out_dir, "ts_sim_dgp_high_min_prpn_0.1_tr_start_24.RData", 
                             sep = ""))
rand_sim_high <- readRDS(paste(out_dir, "rand_sim_dgp_high_min_prpn_0.005_tr_start_24.RData", 
                               sep = ""))
rits_sim_high <- readRDS(paste(out_dir, "rits_sim_dgp_high_min_prpn_0.1_tr_start_24.RData", 
                               sep = ""))

K <- length(unique(rand_sim_high[[1]]$trt))
sim_dat <- readRDS("metadata/sim_dat.RData")

# High SNR
mu_true <- sim_dat$mu_true * 2
contr_true <- mu_true - mu_true[1]
contr_true <- contr_true[setdiff(1:K, sim_dat$placebo_arm)]

summ_rand <- gen_summary_for_table_stoptime(sim = rand_sim_high, K = K, 
                                            ate_ind = ate_ind, contr_true = contr_true, 
                                            need_std = TRUE, need_ipw = FALSE)
summ_ts <- gen_summary_for_table_stoptime(sim = ts_sim_high, K = K,
                                          ate_ind = ate_ind, contr_true = contr_true,
                                          need_ipw = FALSE)
summ_rits <- gen_summary_for_table_stoptime(sim = rits_sim_high, K = K, 
                                            ate_ind = ate_ind, contr_true = contr_true,
                                            need_ipw = FALSE)

cummiscov_rand <- get_cum_mis_cov_stoptime(sim = rand_sim_high, mu_true = mu_true, 
                                           contr_true = contr_true, delay_aipw = 0, 
                                           need_std = TRUE)
cummiscov_ts <- get_cum_mis_cov_stoptime(sim = ts_sim_high, mu_true = mu_true, 
                                         contr_true = contr_true, delay_aipw = 0)
cummiscov_rits <- get_cum_mis_cov_stoptime(sim = ts_sim_high, mu_true = mu_true, 
                                           contr_true = contr_true, delay_aipw = 0)

metric_tab_high <- gen_metric_tab_stoptime(summ_rand = summ_rand, 
                                           summ_ts = summ_ts,
                                           summ_rits = summ_rits, 
                                           cummiscov_rand = cummiscov_rand, 
                                           cummiscov_ts = cummiscov_ts,
                                           cummiscov_rits = cummiscov_rits)
saveRDS(metric_tab_high, "tables/metric_tab_high.RData")
xtable::xtable(metric_tab_high)

# Low SNR
ts_sim_low <- readRDS(paste(out_dir, "ts_sim_dgp_low_min_prpn_0.1_tr_start_24.RData", 
                            sep = ""))
rand_sim_low <- readRDS(paste(out_dir, "rand_sim_dgp_low_min_prpn_0.005_tr_start_24.RData", 
                              sep = ""))
rits_sim_low <- readRDS(paste(out_dir, "rits_sim_dgp_low_min_prpn_0.1_tr_start_24.RData", 
                              sep = ""))

# Low SNR
mu_true <- sim_dat$mu_true
contr_true <- mu_true - mu_true[1]
contr_true <- contr_true[setdiff(1:K, sim_dat$placebo_arm)]

summ_rand <- gen_summary_for_table_stoptime(sim = rand_sim_low, K = K, 
                                            ate_ind = ate_ind, contr_true = contr_true, 
                                            need_std = TRUE, need_ipw = FALSE)
summ_ts <- gen_summary_for_table_stoptime(sim = ts_sim_low, K = K,
                                          ate_ind = ate_ind, contr_true = contr_true,
                                          need_ipw = FALSE)
summ_rits <- gen_summary_for_table_stoptime(sim = rits_sim_low, K = K, 
                                            ate_ind = ate_ind, contr_true = contr_true,
                                            need_ipw = FALSE)

cummiscov_rand <- get_cum_mis_cov_stoptime(sim = rand_sim_low, mu_true = mu_true, 
                                           contr_true = contr_true, delay_aipw = 0, 
                                           need_std = TRUE)
cummiscov_ts <- get_cum_mis_cov_stoptime(sim = ts_sim_low, mu_true = mu_true, 
                                         contr_true = contr_true, delay_aipw = 0)
cummiscov_rits <- get_cum_mis_cov_stoptime(sim = ts_sim_low, mu_true = mu_true, 
                                           contr_true = contr_true, delay_aipw = 0)

metric_tab_low <- gen_metric_tab_stoptime(summ_rand = summ_rand, 
                                          summ_ts = summ_ts,
                                          summ_rits = summ_rits, 
                                          cummiscov_rand = cummiscov_rand, 
                                          cummiscov_ts = cummiscov_ts,
                                          cummiscov_rits = cummiscov_rits)
saveRDS(metric_tab_low, "tables/metric_tab_low.RData")
xtable::xtable(metric_tab_low)

# Null Efficacy
ts_sim_null <- readRDS(paste(out_dir, "ts_sim_dgp_null_min_prpn_0.1_tr_start_24.RData", 
                             sep = ""))
rand_sim_null <- readRDS(paste(out_dir, "rand_sim_dgp_null_min_prpn_0.005_tr_start_24.RData", 
                               sep = ""))
rits_sim_null <- readRDS(paste(out_dir, "rits_sim_dgp_null_min_prpn_0.1_tr_start_24.RData", 
                               sep = ""))

# Null Efficacy
mu_true <- rep(sim_dat$mu_true[1], length(sim_dat$mu_true))
contr_true <- mu_true - mu_true[1]
contr_true <- contr_true[setdiff(1:K, sim_dat$placebo_arm)]

summ_rand <- gen_summary_for_table_stoptime(sim = rand_sim_null, K = K, 
                                            ate_ind = ate_ind, contr_true = contr_true, 
                                            need_std = TRUE, need_ipw = FALSE)
summ_ts <- gen_summary_for_table_stoptime(sim = ts_sim_null, K = K,
                                          ate_ind = ate_ind, contr_true = contr_true,
                                          need_ipw = FALSE)
summ_rits <- gen_summary_for_table_stoptime(sim = rits_sim_null, K = K, 
                                            ate_ind = ate_ind, contr_true = contr_true,
                                            need_ipw = FALSE)

cummiscov_rand <- get_cum_mis_cov_stoptime(sim = rand_sim_null, mu_true = mu_true, 
                                           contr_true = contr_true, delay_aipw = 0, 
                                           need_std = TRUE)
cummiscov_ts <- get_cum_mis_cov_stoptime(sim = ts_sim_null, mu_true = mu_true, 
                                         contr_true = contr_true, delay_aipw = 0)
cummiscov_rits <- get_cum_mis_cov_stoptime(sim = ts_sim_null, mu_true = mu_true, 
                                           contr_true = contr_true, delay_aipw = 0)

metric_tab_null <- gen_metric_tab_stoptime(summ_rand = summ_rand, 
                                           summ_ts = summ_ts,
                                           summ_rits = summ_rits, 
                                           cummiscov_rand = cummiscov_rand, 
                                           cummiscov_ts = cummiscov_ts,
                                           cummiscov_rits = cummiscov_rits)
saveRDS(metric_tab_null, "tables/metric_tab_null.RData")
xtable::xtable(metric_tab_null)

metric_tab_high <- readRDS("tables/metric_tab_high.RData")
metric_tab_low <- readRDS("tables/metric_tab_low.RData")
metric_tab_null <- readRDS("tables/metric_tab_null.RData")
metric_tab <- rbind(metric_tab_high, metric_tab_low, metric_tab_null)
