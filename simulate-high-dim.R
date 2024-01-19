# Perform simulation for one scenario.
library(optparse)
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

# Global parameters
gamma <- 1
sigma <- diag(opt$d)
mu <- rep(0, opt$d)

# Set values for p and k
p <- switch(opt$p,
            low = 2 / opt$n,
            medium = 1 / opt$n,
            high = 1 / (2 * opt$n),
            rlang::abort("Invalid value of p")
)

k <- switch(opt$k,
            large = 0.2 * opt$n,
            medium = 0.1 * opt$n,
            small = 0.05 * opt$n,
            rlang::abort("Invalid value of k")
)

simulate <- function(i) {
  # Generate sample
  data <- mvtnorm::rmvt(opt$n, sigma, 1 / gamma, mu)
  
  # Estimate scatter
  est <- robustbase::covMcd(data, alpha = 0.5)
  
  # Estimate extreme value index and extreme quantile of generating variate
  radius <- sqrt(stats::mahalanobis(data, FALSE, est$cov, inverted = FALSE))
  radius_sort <- sort(radius, decreasing = FALSE)
  
  gamma_est <- mean((log(radius_sort[(opt$n - k):opt$n])
                - log(radius_sort[opt$n - k]))[-1])
  r_hat <- radius_sort[opt$n - k] * (k / (opt$n * p))^gamma_est
  
  # Compute conservative estimate for the error
  compute_error_elliptical(est$cov, r_hat, 1 / gamma, p)
}

# Compute median error from s repetitions and the number of NAs
set.seed(opt$seed)
error <- purrr::map_dbl(1:opt$s, simulate)

num_na <- sum(is.na(error))
error_med <- median(error, na.rm = TRUE)
