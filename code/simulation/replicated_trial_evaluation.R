print(timestamp())
source("code/function/classical_ci.R")
source("code/function/asymp_cs.R")

sim_choice <- readRDS("metadata/sim_choice.RData")
sim_dat <- readRDS("metadata/sim_dat.RData")
dgps <- c("low", "high")
tr_starts <- sim_choice$tr_start
min_prpns <- sim_choice$min_prpns
cases <- expand.grid(dgp = dgps, min_prpn = min_prpns, tr_start = tr_starts)

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
    out_sim <- readRDS(file_name)
    counter <- 0
    if(is.null(out_sim[[1]]$contr)){
      out_sim <- add_asympcs_sim(out_list = out_sim, ate_start = 20, batch = 5, 
                                 placebo_arm = 1, alpha = 0.05, first_peek = 50, 
                                 n_cores = 1)
      counter <- counter + 1
    }
    if(is.null(out_sim[[1]]$contr_standard)){
      out_sim <- add_standard_ci(out = out_sim, ate_start = 20, batch = 5, 
                                 placebo_arm = 1, alpha = 0.05)
      counter <- counter + 1
    }
    if(counter > 0){
      saveRDS(out_sim, file_name)
    } 
  }
  return(NULL)
}, mc.cores = num_cores, cases = cases, sim_choice = sim_choice)
print(timestamp())