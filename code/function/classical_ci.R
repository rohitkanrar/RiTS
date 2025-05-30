standard_ci <- function(rwd_hist, trt_hist, K, placebo_arm = 1, alpha = 0.05){
  # browser()
  if(sum(trt_hist == placebo_arm) == 0){
    out <- matrix(NA, K-1, 3)
  } else{
    n_trt <- as.numeric(table(trt_hist))
    out <- sapply(setdiff(1:K, placebo_arm), function(k){
      if(sum(trt_hist == k) == 0){
        c(NA, NA, NA)
      } else if(sum(trt_hist == placebo_arm) == 1 || sum(trt_hist == k) == 1){
        c(mean(rwd_hist[trt_hist == k]) - mean(rwd_hist[trt_hist == placebo_arm]),
          NA, NA)
      } else{
        tmp <- t.test(rwd_hist[trt_hist == k], rwd_hist[trt_hist == placebo_arm], 
                      var.equal = FALSE)
        c(tmp$estimate[1] - tmp$estimate[2], tmp$conf.int)
      }
    })
    out <- t(out)
  }
  
  dimnames(out) <- list(Arms = setdiff(1:K, placebo_arm), 
                        CI = c("Center", "Lower", "Upper"))
  return(as.matrix(out))
}

add_standard_ci <- function(out, ate_start, batch = 5, placebo_arm = 1, 
                            alpha = 0.05, force_compute = FALSE){
  # browser()
  n_iter <- length(out)
  K <- length(unique(out[[1]]$trt)); N <- length(out[[1]]$trt)
  times <- seq(ate_start, N, batch)
  for(i in 1:n_iter){
    if(is.null(out[[i]]$contr_standard) || force_compute){
      ci_tmp <- array(NA, dim = c(length(times), K-1, 3))
      if(is.null(out[[i]]$reward_benf)){
        reward_benf <- out[[i]]$reward
      } else{
        reward_benf <- out[[i]]$reward_benf
      }
      for(j in 1:length(times)){
        t <- times[j]
        ci_tmp[j, , ] <- standard_ci(rwd_hist = reward_benf[1:t], 
                                     trt_hist = out[[i]]$trt[1:t], K = K,
                                     placebo_arm = placebo_arm, alpha = alpha/(K-1)) # Bonferroni's correction
      }
      dimnames(ci_tmp) <- list(Times = times,
                               Arms = setdiff(1:K, placebo_arm), 
                               CI = c("Center", "Lower", "Upper"))
      out[[i]][["contr_standard"]] <- ci_tmp
    }
  }
  return(out)
}

# rand_sim <- readRDS(paste("output/rand_sim_", 1, ".RData", sep = ""))
# rand_sim1 <- add_standard_ci(rand_sim, times = 50:200)
# saveRDS(rand_sim1, paste("output/rand_sim_", 1, ".RData", sep = ""))
# rm(rand_sim1)
