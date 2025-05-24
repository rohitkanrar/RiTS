source("code/visualization/viz_function.R")
# High SNR
ts_sim_high <- readRDS("output/ts_sim_dgp_high_min_prpn_0.05_tr_start_24.RData")
rand_sim_high <- readRDS("output/rand_sim_dgp_high_min_prpn_0.005_tr_start_24.RData")
rits_sim_high <- readRDS("output/rits_sim_dgp_high_min_prpn_0.05_tr_start_24.RData")
# Low SNR
ts_sim_low <- readRDS("output/ts_sim_dgp_low_min_prpn_0.05_tr_start_24.RData")
rand_sim_low <- readRDS("output/rand_sim_dgp_low_min_prpn_0.005_tr_start_24.RData")
rits_sim_low <- readRDS("output/rits_sim_dgp_low_min_prpn_0.05_tr_start_24.RData")

K <- length(unique(rand_sim_high[[1]]$trt))
sim_choice <- readRDS("metadata/sim_choice.RData")
ate_start <- sim_choice$ate_start
ind <- c(ate_start, seq(30, sim_choice$N, 10))
ate_ind <- sapply(ind, function(i){
  which(as.numeric(dimnames(ts_sim_high[[1]]$contr)[[1]]) == i)
})
sim_dat <- readRDS("metadata/sim_dat.RData")

# High SNR
contr_true <- sim_dat$ite_true_contr * 2

summ_rand <- gen_summary_for_table(sim = rand_sim_high, K = K, 
                                   ate_ind = ate_ind, contr_true = contr_true, 
                                   need_std = TRUE)
summ_ts <- gen_summary_for_table(sim = ts_sim_high, K = K, 
                                   ate_ind = ate_ind, contr_true = contr_true)
summ_rits <- gen_summary_for_table(sim = rits_sim_high, K = K, 
                                   ate_ind = ate_ind, contr_true = contr_true)

est_err_tab_high <- gen_bias_rmse_tab(summ_rand = summ_rand, 
                                      summ_ts = summ_ts,
                                      summ_rits = summ_rits, 
                                      ate_ind = ate_ind, ind = ind)
xtable::xtable(est_err_tab_high)

# Low SNR
contr_true <- sim_dat$ite_true_contr

summ_rand <- gen_summary_for_table(sim = rand_sim_low, K = K, 
                                   ate_ind = ate_ind, contr_true = contr_true, 
                                   need_std = TRUE)
summ_ts <- gen_summary_for_table(sim = ts_sim_low, K = K, 
                                 ate_ind = ate_ind, contr_true = contr_true)
summ_rits <- gen_summary_for_table(sim = rits_sim_low, K = K, 
                                   ate_ind = ate_ind, contr_true = contr_true)

est_err_tab_low <- gen_bias_rmse_tab(summ_rand = summ_rand, 
                                      summ_ts = summ_ts,
                                      summ_rits = summ_rits, 
                                      ate_ind = ate_ind, ind = ind)
xtable::xtable(est_err_tab_low)
