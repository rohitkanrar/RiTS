source("code/visualization/functions.R")
source("code/function/misc.R")
library(reshape2); library(ggplot2); library(patchwork)
sim_choice <- readRDS("metadata/sim_choice.RData")
sim_dat <- readRDS("metadata/sim_dat.RData")
dgps <- c("low", "high")
tr_starts <- sim_choice$tr_start
min_prpns <- sim_choice$min_prpns
cases <- expand.grid(dgp = dgps, min_prpn = min_prpns, tr_start = tr_starts)
out_dir <- "output/varyingX/"
alpha <- 0.05; ate_start <- sim_choice$ate_start; K <- 4

cum_miscov_plot <- vector("list", 2*nrow(cases))
j <- 1
for(i in 1:nrow(cases)){
  dgp <- cases[i, "dgp"]; min_prpn <- cases[i, "min_prpn"]
  tr_start <- cases[i, "tr_start"]
  case_str <- paste("dgp", dgp, "min_prpn", min_prpn, "tr_start", tr_start, 
                    sep = "_")
  # mods <- c("rand_sim", "ts_sim", "rits_sim", "rand_mis_sim", "ts_mis_sim", 
  #           "rits_mis_sim")
  mods <- c("ts_sim", "rits_sim", "ts_mis_sim", "rits_mis_sim")
  for(mod in mods){
    if(dgp == "high"){
      mu_true <- sim_dat$mu_true * 2
    } else{
      mu_true <- sim_dat$mu_true
    }
    contr_true <- mu_true - mu_true[1]
    contr_true <- contr_true[setdiff(1:K, sim_dat$placebo_arm)]
    ts_sim <- readRDS(paste(out_dir, "ts_sim_", case_str, ".RData", sep = ""))
    rits_sim <- readRDS(paste(out_dir, "rits_sim_", case_str, ".RData", sep = ""))
    tmp1 <- gen_cum_miscov_plot(out = ts_sim, ate_true = mu_true, 
                                              contr_true = contr_true, 
                                              alpha = alpha, ate_start = ate_start,
                                              titl = "TS")
    tmp2 <- gen_cum_miscov_plot(out = rits_sim, ate_true = mu_true, 
                                                contr_true = contr_true, 
                                                alpha = alpha, ate_start = ate_start,
                                                titl = "RiTS")
    tmp <- tmp1 + tmp2 + plot_layout(ncol = 2, guides = "collect")
    flname <- paste("plot/cum_miscov_all/", case_str, ".jpg", sep = "")
    ggsave(flname, tmp, height = 4, width = 12, units = "in")
  }
}