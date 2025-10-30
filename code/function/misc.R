sigmoid <- function(x) 1/(1 + exp(-x))

# y = a - b[(x+c)^2 + (x-c)^2]
get_beta <- function(a, b, c) c(a-2*b*c^2, 0, -2*b)
# y = a - b x^2
get_gamma <- function(a, b) c(a, 0, -b)

get_X <- function(N){
  x <- rnorm(N)
  cbind(1, x, x^2)
}

get_true_avg_rwd <- function(X, trt, beta_true){
  # browser()
  N <- nrow(X)
  K <- length(unique(trt))
  mu <- matrix(0, N, K)
  for(n in 1:N){
    mu[n, ] <- as.numeric(crossprod(X[n, ], beta_true))
  }
  mu
}

apply_floor <- function(a, amin) {
  # browser()
  new <- pmax(a, amin)
  total_slack <- sum(new) - 1
  individual_slack <- new - amin
  c <- total_slack / sum(individual_slack)
  return(new - c * individual_slack)
}

expand_standard_array <- function(target_array, standard_array){
  target_times_chr <- dimnames(target_array)[[1]]
  standard_times_chr <- dimnames(standard_array)[[1]]
  target_times_num <- as.numeric(target_times_chr)
  standard_times_num <- as.numeric(standard_times_chr)
  source_indices <- findInterval(target_times_num, standard_times_num)
  source_indices[source_indices == 0] <- 1
  new_array <- standard_array[source_indices, , , drop = FALSE]
  dimnames(new_array) <- dimnames(target_array)
  return(new_array)
}

get_cum_mis_cov <- function(sim, mu_true, contr_true, delay_aipw = 0, 
                            delay_ipw = 0, need_std = FALSE, need_ipw = FALSE){
  n_iter <- length(sim)
  total_peek <- dim(sim[[1]]$ate)[1]
  no_of_peek <- dim(sim[[1]]$ate)[1] - delay_aipw
  cum_mis_cov_ate <- matrix(0, no_of_peek, dim(sim[[1]]$ate)[2])
  cum_mis_cov_contr <- matrix(0, no_of_peek, dim(sim[[1]]$contr)[2])
  if(need_std){
    cum_mis_cov_contr_std <- matrix(0, no_of_peek, dim(sim[[1]]$contr)[2])
  }
  if(need_ipw){
    no_of_peek_ipw <- dim(sim[[1]]$ate)[1] - delay_ipw
    cum_mis_cov_contr_ipw <- matrix(0, no_of_peek_ipw, dim(sim[[1]]$contr)[2])
  }
  for(iter in 1:n_iter){
    sim[[iter]][["contr_standard"]] <- expand_standard_array(target_array = sim[[iter]][["contr"]],
                                                             standard_array = sim[[iter]][["contr_standard"]])
    for(k in 1:K){
      cum_mis_cov_ate[, k] <- cum_mis_cov_ate[, k] + 
        cummax(sim[[iter]]$ate[(delay_aipw+1):total_peek, k,  2] > mu_true[k] | 
                 sim[[iter]]$ate[(delay_aipw+1):total_peek, k, 3] < mu_true[k])
      if(k > 1){
        cum_mis_cov_contr[, k-1] <-   cum_mis_cov_contr[, k-1] + 
          cummax(sim[[iter]]$contr[(delay_aipw+1):total_peek, k-1,  2] > contr_true[k-1] | 
                   sim[[iter]]$contr[(delay_aipw+1):total_peek, k-1, 3] < contr_true[k-1])
        if(need_std){
          cum_mis_cov_contr_std[, k-1] <-   cum_mis_cov_contr_std[, k-1] + 
            cummax(sim[[iter]]$contr_standard[(delay_aipw+1):total_peek, k-1,  2] > contr_true[k-1] | 
                     sim[[iter]]$contr_standard[(delay_aipw+1):total_peek, k-1, 3] < contr_true[k-1])
        }
        if(need_ipw){
          cum_mis_cov_contr_ipw[, k-1] <-   cum_mis_cov_contr_ipw[, k-1] + 
            cummax(sim[[iter]]$contr_ipw[(delay_ipw+1):total_peek, k-1,  2] > contr_true[k-1] | 
                     sim[[iter]]$contr_ipw[(delay_ipw+1):total_peek, k-1, 3] < contr_true[k-1])
        }
      } 
    }
  }
  out_list <- list(ate = cum_mis_cov_ate, contr = cum_mis_cov_contr)
  if(need_std){
    out_list[["contr_std"]] <- cum_mis_cov_contr_std
  } 
  if(need_ipw){
    out_list[["contr_ipw"]] <- cum_mis_cov_contr_ipw
  }
  return(out_list)
}

zero_in_intv <- function(intv, zero = 0.1){
  a <- intv[1]; b <- intv[2]
  if(is.na(a) || is.na(b)){
    return(FALSE)
  } else if(a < zero && b > zero){
    return(TRUE)
  } else{
    return(FALSE)
  }
}

zero_where_intv <- function(intv, zero = 0.1){
  a <- intv[1]; b <- intv[2]
  if(is.na(a) || is.na(b)){
    return("NOIDEA")
  } else if(b < zero){
    return("ABOVE")
  } else if(a > zero){
    return("BELOW")
  } else if(a < zero && b > zero){
    return("INSIDE")
  } else{
    return("NOIDEA")
  }
}

stop_trial <- function(intvs, K, zero = 0.1){
  pos <- sapply(1:(K-1), function(k){
    zero_where_intv(intvs[k, ], zero = zero)
  })
  if("BELOW" %in% pos){
    return(1)
  } else if(sum(pos == "ABOVE") == (K-1)){
    return(1)
  } else{
    return(0)
  }
}
