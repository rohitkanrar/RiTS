source("code/visualization/viz_function.R")
sim_choice <- readRDS("metadata/sim_choice.RData")
sim_dat <- readRDS("metadata/sim_dat.RData")
dgps <- c("low", "high")
tr_starts <- sim_choice$tr_start
min_prpns <- sim_choice$min_prpns
cases <- expand.grid(dgp = dgps, min_prpn = min_prpns, tr_start = tr_starts)
n_iter <- 5000; out_dir <- "output/"; K <- sim_choice$K; N <- sim_choice$N
ind <- round(
  c(seq(50, sim_choice$N*0.625, 10), sim_choice$N*(3/4), sim_choice$N/2)
)
j <- 0

# cumulative mis-coverage and estimation error tables
cum_miscov_all_tab <- matrix(NA, nrow(cases), 4*(K-1))
cum_miscov_all_mis_tab <- matrix(NA, nrow(cases), 4*(K-1))
est_err_all_tab <- vector(mode = "list", length = nrow(cases))
est_err_all_mis_tab <- vector(mode = "list", length = nrow(cases))

for(i in 1:nrow(cases)){
  print(i)
  dgp <- cases[i, "dgp"]; min_prpn <- cases[i, "min_prpn"]
  tr_start <- cases[i, "tr_start"]
  case_str <- paste("dgp", dgp, "min_prpn", min_prpn, "tr_start", tr_start, 
                    sep = "_")
  case_str0 <- paste("dgp", dgp, "min_prpn", 0.005, "tr_start", sim_choice$ate_start, 
                     sep = "_")
  rand_file_name <- paste(out_dir, "rand_sim_", case_str0, ".RData", sep = "")
  rand_mis_file_name <- paste(out_dir, "rand_mis_sim_", case_str0, ".RData", sep = "")
  ts_file_name <- paste(out_dir, "ts_sim_", case_str, ".RData", sep = "")
  ts_mis_file_name <- paste(out_dir, "ts_mis_sim_", case_str, ".RData", sep = "")
  rits_file_name <- paste(out_dir, "rits_sim_", case_str, ".RData", sep = "")
  rits_mis_file_name <- paste(out_dir, "rits_mis_sim_", case_str, ".RData", sep = "")
  
  if(dgp == "high"){
    mu_true <- sim_dat$mu_true * 2
  } else{
    mu_true <- sim_dat$mu_true
  }
  contr_true <- mu_true - mu_true[1]
  contr_true <- contr_true[setdiff(1:K, sim_dat$placebo_arm)]
  ## Correctly specified
  # cumulative mis-coverage
  rand_sim <- readRDS(rand_file_name); rand_sim <- rand_sim[1:n_iter]
  if(j == 0){
    ate_ind <- sapply(ind, function(i){
      which(as.numeric(dimnames(rand_sim[[1]]$contr)[[1]]) == i)
    })
    j <- j + 1
  }
  ts_sim <- readRDS(ts_file_name); ts_sim <- ts_sim[1:n_iter]
  rits_sim <- readRDS(rits_file_name); rits_sim <- rits_sim[1:n_iter]
  cum_miscov_all_tab[i, ] <- gen_cum_miscov_for_tab(out_rand = rand_sim, 
                                                    out_ts = ts_sim, 
                                                    out_rits = rits_sim, 
                                                    mu_true = mu_true, 
                                                    contr_true = contr_true)
  # estimation error (bias and rmse)
  summ_rand <- gen_summary_for_table(sim = rand_sim, K = K, 
                                     ate_ind = ate_ind, contr_true = contr_true, 
                                     need_std = TRUE)
  summ_ts <- gen_summary_for_table(sim = ts_sim, K = K, 
                                   ate_ind = ate_ind, contr_true = contr_true)
  summ_rits <- gen_summary_for_table(sim = rits_sim, K = K, 
                                     ate_ind = ate_ind, contr_true = contr_true)
  
  est_err_tab <- gen_bias_rmse_tab(summ_rand = summ_rand, 
                                       summ_ts = summ_ts,
                                       summ_rits = summ_rits, 
                                       ate_ind = ate_ind, ind = ind)
  est_err_all_tab[[i]] <- list(tab = est_err_tab, case_str = case_str)
  ## Mis-specified
  # cumulative mis-coverage
  rand_sim <- readRDS(rand_mis_file_name); ts_sim <- readRDS(ts_mis_file_name); 
  rits_sim <- readRDS(rits_mis_file_name)
  cum_miscov_all_mis_tab[i, ] <- gen_cum_miscov_for_tab(out_rand = rand_sim, 
                                                    out_ts = ts_sim, 
                                                    out_rits = rits_sim, 
                                                    mu_true = mu_true, 
                                                    contr_true = contr_true)
  # estimation error (bias and rmse)
  summ_rand <- gen_summary_for_table(sim = rand_sim, K = K, 
                                     ate_ind = ate_ind, contr_true = contr_true, 
                                     need_std = TRUE)
  summ_ts <- gen_summary_for_table(sim = ts_sim, K = K, 
                                   ate_ind = ate_ind, contr_true = contr_true)
  summ_rits <- gen_summary_for_table(sim = rits_sim, K = K, 
                                     ate_ind = ate_ind, contr_true = contr_true)
  
  est_err_tab <- gen_bias_rmse_tab(summ_rand = summ_rand, 
                                   summ_ts = summ_ts,
                                   summ_rits = summ_rits, 
                                   ate_ind = ate_ind, ind = ind)
  est_err_all_mis_tab[[i]] <- list(tab = est_err_tab, case_str = case_str)
}

saveRDS(cum_miscov_all_tab, "tables/cum_miscov_all_tab.RData")
saveRDS(cum_miscov_all_mis_tab, "tables/cum_miscov_all_mis_tab.RData")
saveRDS(est_err_all_tab, "tables/est_err_all_tab.RData")
saveRDS(est_err_all_mis_tab, "tables/est_err_all_mis_tab.RData")