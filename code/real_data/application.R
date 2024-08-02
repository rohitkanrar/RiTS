real_dat <- readRDS("data/clean_data.RData")
countfact_eff <- readRDS("output/countfact_eff.RData")
countfact_safe <- readRDS("output/countfact_safe.RData")

X_real <- as.matrix(real_dat[c("baseline_salt", "baseline_lymph", "age", "bmi",
                               "albumin", "lymph_marker")])
X_real <- cbind(1, X_real, (as.numeric(real_dat$severity)-1))

source("code/function/main_function.R")
source("code/real_data/application_function.R")
source("code/function/misc.R")
source("code/function/awaipw_functions.R")
set.seed(2024)

K <- length(unique(real_dat[["arm"]]))
v <- 10
ts_real <- do_ts_batch_real(X_real, countfact_eff, countfact_safe, K, v = v)
rits_real <- do_rits_batch_real(X_real, countfact_eff, countfact_safe, K, 1, v = v)
rand_real <- get_metrics_rand(countfact_eff, countfact_safe, real_dat[["arm"]])

saveRDS(ts_real, "output/ts_real.RData")
saveRDS(rits_real, "output/rits_real.RData")
saveRDS(rand_real, "output/rand_real.RData")

rm(list = ls())
