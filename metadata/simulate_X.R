set.seed(479612)
sim_choice <- readRDS("metadata/sim_choice.RData")
sim_dat <- readRDS("metadata/sim_dat.RData")
N <- sim_choice$N;
x <- rnorm(N)
X <- cbind(1, x)
X_true <- cbind(1, x, x^2)
sim_dat[["X"]] <- X
sim_dat[["X_true"]] <- X_true
saveRDS(sim_dat, "metadata/sim_dat.RData")

if(FALSE){
  ite_true <- sim_dat$X_true %*% sim_dat$beta_true[, , 1]
  ite_true_mean <- sapply(1:4, function(j) cumsum(ite_true[, j])/(1:100))
  
  plot(ite_true_mean[, 1], type = "l", ylim = c(0.9, 1.45))
  lines(ite_true_mean[, 2]); lines(ite_true_mean[, 3]); lines(ite_true_mean[, 4])
  abline(h = sim_dat$mu_true[1], col = "red")
  abline(h = sim_dat$mu_true[2], col = "red")
  abline(h = sim_dat$mu_true[3], col = "red")
  abline(h = sim_dat$mu_true[4], col = "red")
}

