# source("code/function/classical_ci.R")
source("code/function/asymp_cs.R")
require(parallel)
rand_files <- list.files("output/varyingX/", pattern = "^rand_")
mclapply(1:length(rand_files), function(i){
  flname <- paste("output/varyingX/", rand_files[i], sep = "")
  batch <- 1
  out_sim <- readRDS(flname)
  # out_sim <- add_standard_ci(out = out_sim, ate_start = 20, batch = batch, 
  #                            placebo_arm = 1, alpha = 0.05, force_compute = TRUE)
  out_sim <- add_asympcs_sim(out_list = out_sim, ate_start = 20, batch = batch, 
                             placebo_arm = 1, alpha = 0.05, first_peek = 50, 
                             n_cores = 1, force_compute = TRUE)
  saveRDS(object = out_sim, file = flname)
}, mc.cores = 9)
