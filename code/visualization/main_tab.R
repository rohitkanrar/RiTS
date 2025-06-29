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
ind <- round(
  c(seq(50, sim_choice$N*0.625, 10), sim_choice$N*(3/4), sim_choice$N)
)

ate_ind <- sapply(ind, function(i){
  which(as.numeric(dimnames(ts_sim_high[[1]]$contr)[[1]]) == i)
})
sim_dat <- readRDS("metadata/sim_dat.RData")

# High SNR
mu_true <- sim_dat$mu_true * 2
contr_true <- mu_true - mu_true[1]
contr_true <- contr_true[setdiff(1:K, sim_dat$placebo_arm)]

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
saveRDS(est_err_tab_high, "tables/est_err_tab_high.RData")
xtable::xtable(est_err_tab_high)

# Low SNR
mu_true <- sim_dat$mu_true
contr_true <- mu_true - mu_true[1]
contr_true <- contr_true[setdiff(1:K, sim_dat$placebo_arm)]

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
saveRDS(est_err_tab_low, "tables/est_err_tab_low.RData")
xtable::xtable(est_err_tab_low)

# est_err_tab_high <- readRDS("tables/est_err_tab_high.RData")
# est_err_tab_low <- readRDS("tables/est_err_tab_low.RData")
# est_err_tab <- rbind(est_err_tab_high, est_err_tab_low)