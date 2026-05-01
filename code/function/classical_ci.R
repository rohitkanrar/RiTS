if(FALSE){
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
}



standard_ci <- function(rwd_hist, trt_hist, K, placebo_arm = 1, c_k){
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
                      var.equal = TRUE)
        est <- tmp$estimate[1] - tmp$estimate[2]
        low_ci <- est - c_k * tmp$stderr
        up_ci <- est + c_k * tmp$stderr
        c(est, low_ci, up_ci)
      }
    })
    out <- t(out)
  }
  
  dimnames(out) <- list(Arms = setdiff(1:K, placebo_arm), 
                        CI = c("Center", "Lower", "Upper"))
  return(as.matrix(out))
}

add_standard_ci <- function(out, ate_start, n_looks = 30, 
                            placebo_arm = 1, force_compute = FALSE, c_ks = NULL){
  # browser()
  n_iter <- length(out)
  K <- length(unique(out[[1]]$trt)); N <- length(out[[1]]$trt)
  times <- floor(seq(ate_start, N, length.out = n_looks))
  
  # If c_ks is not supplied, load the pre-calculated group-sequential boundaries
  if (is.null(c_ks)) {
    if (file.exists("metadata/ck.RData")) {
      c_ks <- readRDS("metadata/ck.RData")
    } else {
      stop(
        "c_ks was not supplied and metadata/ck.RData was not found. ",
        "Either create metadata/ck.RData first or pass c_ks explicitly."
      )
    }
  }
  
  # If one critical value is supplied, repeat it for all looks
  if (length(c_ks) == 1) {
    c_ks <- rep(c_ks, length(times))
  }
  
  # Check that c_ks has the correct length
  if (length(c_ks) != length(times)) {
    stop(
      "Length of c_ks must be either 1 or equal to the number of looks. ",
      "length(c_ks) = ", length(c_ks),
      ", length(times) = ", length(times), "."
    )
  }
  
  # converting bounds for known variance to unknown variance
  c_ks_tbounds <- qt(1 - pnorm(c_ks), times - 1, lower.tail = FALSE)
  if (any(is.infinite(c_ks_tbounds))) {
    c_ks_tbounds[is.infinite(c_ks_tbounds)] <- max(c_ks_tbounds[is.finite(c_ks_tbounds)])
  }
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
                                     placebo_arm = placebo_arm, 
                                     c_k = c_ks_tbounds[j]) # O'brien & Fleming ESF with variance unknown
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
