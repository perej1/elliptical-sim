# Perform simulation for one scenario.
source("functions.R")

option_list <- list(
  make_option("--s", type = "integer", default = 10,
              help = "Number of repetitions for a scenario"),
  make_option("--d", type = "integer", default = 2,
              help = "Dimensions"),
  make_option("--n", type = "integer", default = 5000,
              help = "Sample size"),
  make_option("--p", type = "character", default = "high",
              help = "Probability mass outside quantile region"),
  make_option("--k", type = "character", default = "large",
              help = "Sample size of the tail"),
  make_option("--seed", type = "integer", default = 278,
              help = "Set seed for sampling")
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

simulate <- function(s, d, n, p, k, seed) {
  # Global parameters
  gamma <- 1
  sigma <- diag(d)
  mu <- rep(0, d)
  
  # Set values for p and k
  p <- switch(p,
              low = 2 / n,
              medium = 1 / n,
              high = 1 / (2 * n),
              rlang::abort("Invalid value of p")
  )
  
  k <- switch(opt$k,
              large = 0.2 * n,
              medium = 0.1 * n,
              small = 0.05 * n,
              rlang::abort("Invalid value of k")
  )
  
  simulate_round <- function(i) {
    # Generate sample
    data <- mvtnorm::rmvt(n, sigma, 1 / gamma, mu)
    
    # Estimate scatter
    est <- robustbase::covMcd(data, alpha = 0.5)
    
    # Estimate extreme value index and extreme quantile of generating variate
    radius <- sqrt(stats::mahalanobis(data, FALSE, est$cov, inverted = FALSE))
    radius_sort <- sort(radius, decreasing = FALSE)
    
    gamma_est <- mean((log(radius_sort[(n - k):n])
                       - log(radius_sort[n - k]))[-1])
    r_hat <- radius_sort[n - k] * (k / (n * p))^gamma_est
    
    # Compute conservative estimate for the error
    compute_error_elliptical(est$cov, r_hat, 1 / gamma, p)
  }
  
  # Compute median error from s repetitions and the number of NAs
  set.seed(seed)
  error <- purrr::map_dbl(1:s, simulate_round)
  
  list(error_med = median(error, na.rm = TRUE), num_na = sum(is.na(error)))
}

s <- 10
d <- 2
n <- 5000
p <- "high"
k <- "large"
seed <- 278

simulate(s, d, n, p, k, seed)
