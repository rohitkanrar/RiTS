print(timestamp())
source("code/function/classical_ci.R")
source("code/function/asymp_cs.R")

sim_choice <- readRDS("metadata/sim_choice.RData")
dgps <- c("low", "high", "null")
tr_starts <- sim_choice$tr_start
min_prpns <- sim_choice$min_prpns
cases <- expand.grid(dgp = dgps, min_prpn = min_prpns, tr_start = tr_starts)
batch <- 1

library(parallel)
num_cores <- 16
results <- mclapply(1:nrow(cases), function(i, cases, sim_choice){
  # browser()
  out_dir <- "output/"
  dgp <- cases[i, "dgp"]; min_prpn <- cases[i, "min_prpn"]
  tr_start <- cases[i, "tr_start"]
  case_str <- paste("dgp", dgp, "min_prpn", min_prpn, "tr_start", tr_start, 
                    sep = "_")
  mods <- c("rand_sim", "ts_sim", "rits_sim", "rand_mis_sim", "ts_mis_sim", 
            "rits_mis_sim")
  for(mod in mods){
    file_name <- paste(out_dir, mod, "_", case_str, ".RData", sep = "")
    if(file.exists(file_name)){
      out_sim <- readRDS(file_name)
      # out_sim <- add_asympcs_sim(out_list = out_sim, ate_start = 24, batch = batch,
      #                            placebo_arm = 1, alpha = 0.05, first_peek = 100,
      #                            n_cores = 1, force_compute = TRUE,
      #                            learner = "main_ridge")
      # out_sim <- add_asympcs_sim(out_list = out_sim, ate_start = 24, batch = batch,
      #                            placebo_arm = 1, alpha = 0.05, first_peek = 50,
      #                            n_cores = 1, force_compute = TRUE,
      #                            learner = NULL)
      out_sim <- add_standard_ci(out = out_sim, ate_start = 24, n_looks = 30, 
                                 placebo_arm = 1, force_compute = TRUE)
      saveRDS(out_sim, file_name) 
    } else{
      next
    }
  }
  return(NULL)
}, mc.cores = num_cores, cases = cases, sim_choice = sim_choice)
print(timestamp())