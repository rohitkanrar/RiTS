sigmoid <- function(x) 1/(1 + exp(-x))

# y = a - b[(x+c)^2 + (x-c)^2]
get_beta <- function(a, b, c) c(a-2*b*c^2, 0, -2*b)
# y = a - b x^2
get_gamma <- function(a, b) c(a, 0, -b)

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
  new <- pmax(a, amin)
  total_slack <- sum(new) - 1
  individual_slack <- new - amin
  c <- total_slack / sum(individual_slack)
  return(new - c * individual_slack)
}

get_cum_mis_cov <- function(sim, mu_true, contr_true, delay = 0){
  total_peek <- dim(sim[[1]]$ate)[1]
  no_of_peek <- dim(sim[[1]]$ate)[1] - delay
  cum_mis_cov_ate <- matrix(0, no_of_peek, dim(sim[[1]]$ate)[2])
  cum_mis_cov_contr <- matrix(0, no_of_peek, dim(sim[[1]]$contr)[2])
  for(iter in 1:n_iter){
    for(k in 1:K){
      cum_mis_cov_ate[, k] <- cum_mis_cov_ate[, k] + 
        cummax(sim[[iter]]$ate[(delay+1):total_peek, k,  2] > mu_true[k] | 
                 sim[[iter]]$ate[(delay+1):total_peek, k, 3] < mu_true[k])
      if(k > 1){
        cum_mis_cov_contr[, k-1] <-   cum_mis_cov_contr[, k-1] + 
          cummax(sim[[iter]]$contr[(delay+1):total_peek, k-1,  2] > contr_true[k-1] | 
                   sim[[iter]]$contr[(delay+1):total_peek, k-1, 3] < contr_true[k-1])
      } 
    }
  }
  list(cum_mis_cov_ate, cum_mis_cov_contr)
}
