source("code/function/asymp_cs.R")
source("code/function/misc.R")
# ts_sim <- readRDS("output/ts_sim.RData")
rits_sim <- readRDS("output/rits_sim.RData")
sim_dat <- readRDS("output/sim_dat.RData")
repl <- length(rits_sim)
K <- length(sim_dat$mu_true)
N <- nrow(sim_dat$X_true)
m <- 50
beta_true <- sim_dat$beta_true
ts_sim_asympcs <- array(0, c(repl, K-1, N-m+1, 2))
coverage <- array(0, c(repl, N-m+1, K-1))
coverage_static <- array(0, c(repl, N-m+1, K-1))
contr_true <- sim_dat$mu_true - sim_dat$mu_true[1]

for(r in 1:repl){
  print(r)
  # out <- ts_sim[[r]]
  out <- rits_sim[[r]]
  # tmp <- get_asympcs(out$trt, out$reward, out$log_dat$prpns_mat,
  #                    out$log_dat$context, 1, times = m:N,
  #                    first_peek = N/2)
  # tmp <- get_asympcs(out$trt, out$reward_benf, out$log_dat$prpns_mat,
  #                    out$log_dat$context, 1, times = m:N,
  #                    first_peek = N/2)
  true_mu <- get_true_avg_rwd(sim_dat$X_true, out$trt, sim_dat$beta_true[, , 1])
  
  png(paste("plot/asympcs_coverage_rits/", r, ".png", sep = ""),
      width = 800, height = 600)
  par(mfrow = c(1, 3))
  for(k in 2:K){
    # tmp1 <- as.matrix(tmp[[k-1]])
    tmp1 <- ts_sim_asympcs[r, k-1, , ]
    # ts_sim_asympcs[r, k-1, , ] <- tmp1
    true_mu_running <- cumsum(true_mu[, k] - true_mu[, 1]) / (1:N)
    true_mu_running <- true_mu_running[m:N]
    # coverage[r, , k-1] <- cummax(tmp1[, 1] > true_mu_running | 
    #                                tmp1[, 2] < true_mu_running)
    coverage_static[r, , k-1] <- cummax(tmp1[, 1] > contr_true[k] |
                                   tmp1[, 2] < contr_true[k])
    
    plot(ts_sim_asympcs[r, k-1, , 1], ylim = c(-2, 3), type = "l")
    lines(ts_sim_asympcs[r, k-1, , 2])
    lines(true_mu_running, col = "red")
    abline(h = contr_true[k], col = "blue")
  }
  par(mfrow = c(1, 1))
  dev.off()
}
saveRDS(ts_sim_asympcs, "output/ts_sim_asympcs_rits.RData")
saveRDS(coverage, "output/coverage_rits.RData")
saveRDS(coverage_static, "output/coverage_static_rits.RData")

png(paste("plot/asympcs_cum_coverage", r, ".png", sep = ""),
    width = 900, height = 300)
par(mfrow = c(1, 3))
for(k in 2:K){
  plot(m:N, apply(coverage[, , k-1], 2, mean), xlab = "Patient", 
       ylim = c(0, 0.06), ylab = "Cumulative Miscoverage", 
       main = paste("Arm", k, "- Arm 1"), type = "l", col = "red")
  lines(m:N, apply(coverage_static[, , k-1], 2, mean), col = "blue")
  abline(h = 0.05)
}
par(mfrow = c(1, 1))
dev.off()

png(paste("plot/asympcs_cum_coverage_static", r, ".png", sep = ""),
    width = 900, height = 300)
par(mfrow = c(1, 3))
for(k in 2:K){
  plot(m:N, apply(coverage_static[, , k-1], 2, mean), xlab = "Patient", 
       ylim = c(0, 0.06), ylab = "Cumulative Miscoverage", 
       main = paste("Arm", k, "- Arm 1"), type = "l")
  abline(h = 0.05)
}
par(mfrow = c(1, 1))
dev.off()





