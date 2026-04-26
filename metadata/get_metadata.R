sim_choice <- list(
  N = 200,
  K = 4,
  d = 3,
  tr_start = c(24, 40, 56),
  ate_start = 24,
  reward_sig = 1,
  alpha = 0.05,
  design = "clip",
  min_prpns = c(0.005, 0.01, 0.05, 0.1)
)

beta_true <- read.csv("metadata/beta_true.csv")
beta_true_arr <- array(NA, dim = c(3, 4, 2))
beta_true_arr[, , 1] <- as.matrix(beta_true[, 2:5])
beta_true_arr[, , 2] <- as.matrix(beta_true[, 6:9])

X_true <- read.csv("metadata/X_true.csv")
X <- X_true[, -3]
sim_dat <- list(
  X = X, 
  X_true = X_true,
  mu_true = c(0.988, 1.1, 1.225, 1.35),
  nu_true = c(1, 0.995, 0.95, 0.7),
  util_true = c(0.994, 1.047, 1.087, 1.025),
  beta_true = beta_true_arr,
  placebo_arm = 1, 
  weight = 1
)

saveRDS(sim_dat, "metadata/sim_dat.RData")
saveRDS(sim_choice, "metadata/sim_choice.RData")
rm(list = ls())


library(gsDesign)
k <- 30
info.frac <- (1:k) / k
design <- gsDesign(
  k = k,
  n.I = info.frac,
  test.type = 1,
  alpha = 0.05/(3 * 2), # 3 for number of active doses, 2 for two-sided CI
  sfu = sfLDOF
)
c_k <- design$upper$bound
saveRDS(c_k, "metadata/ck.RData")
rm(list = ls())