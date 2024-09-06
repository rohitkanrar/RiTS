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
