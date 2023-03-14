library(optparse)
library(mvtnorm)
suppressPackageStartupMessages(library(dplyr))


get_ball_mesh <- function(m) {
  w <- seq(0, 2 * pi, length.out = m)
  cbind(cos(w), sin(w))
}


sqrtmat <- function(sigma) {
  eigenval <- eigen(sigma)$values
  if (any(eigenval <= 0) || any(sigma != t(sigma))) {
    rlang::abort("`sigma` must be a symmetric positive definite matrix.")
  }
  eigenvec <- eigen(sigma)$vectors
  eigenvec %*% diag(eigenval^0.5) %*% t(eigenvec)
}


tdist_extreme_region <- function(sigma, gamma, p, m) {
  w <- get_ball_mesh(m)
  lambda <- sqrtmat(d * stats::qf(1 - p, d, 1 / gamma) * sigma)
  w %*% t(lambda)
}


elliptical_extreme_qregion <- function(data, sigma_est, p, k, m) {
  n <- nrow(data)
  w <- get_ball_mesh(m)
  
  # Approximate generating variate
  radius <- sqrt(stats::mahalanobis(data, FALSE, sigma_est, inverted = FALSE))
  radius_sort <- sort(radius, decreasing = FALSE)
  
  # Estimate extreme value index
  gamma_est <- mean((log(radius_sort[(n - k):n]) - log(radius_sort[n - k]))[-1])
  
  # Estimate extreme quantile of generating variate
  r_hat <- radius_sort[n - k] * (k / (n * p))^gamma_est
  
  # Estimate extreme quantile region
  lambda <- sqrtmat(r_hat^2 * sigma_est)
  w %*% t(lambda)
}


depth_extreme_qregion <- function(data, p, k, m) {
  n <- nrow(data)
  w <- get_ball_mesh(m)
  
  radius <- apply(data, 1, norm, type = "2")
  radius_sort <- sort(radius, decreasing = FALSE)
  
  gamma_est <- mean((log(radius_sort[(n - k):n]) - log(radius_sort[n - k]))[-1])
  u_est <- radius_sort[n - k]
  nu_hat_unit_halfspace <- 1 / k * apply(w %*% t(data) >= u_est, 1, sum)
  hd_w_nu_hat_star <- sweep(pmax(w %*% t(w), 0)^(-1 / gamma_est), 2,
                            nu_hat_unit_halfspace, "*")
  hd_w_nu_hat_star <- apply(hd_w_nu_hat_star, 1, min, na.rm = TRUE)
  
  data_unit <- sweep(data, 1, radius, "/")
  x_approx_w <- apply(data_unit %*% t(w), 1, which.max)
  nu_s_hat <- 1 / k * sum(radius / u_est >= hd_w_nu_hat_star[x_approx_w]^gamma_est)
  
  r <- u_est * (k * nu_s_hat / (n * p))^gamma_est * hd_w_nu_hat_star^gamma_est
  r * w
}


option_list <- list(
  make_option("--type", type = "character", default = "cauchy",
              help = "Distribution type"),
  make_option("--n", type = "integer", default = 500,
              help = "Sample size"),
  make_option("--seed", type = "integer", default = 204,
              help = "Set seed for sampling"),
  make_option("--p", type = "character", default = "low",
              help = "Probability mass outside quantile region"),
  make_option("--k", type = "character", default = "large",
              help = "Sample size of the tail")
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Global parameters
s <- 100
d <- 2
m <- 1000

# Scatter and location
if (opt$type == "cauchyAff") {
  sigma <- matrix(c(11, 10.5, 10.5, 11.25), byrow = TRUE, ncol = d)
} else {
  sigma <- diag(d)
}
mu <- rep(0, d)

# Set values for p, k and extreme value index
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

gamma <- switch(opt$type,
                cauchy = 1,
                cauchyAff = 1,
                cauchy3d = 1,
                tdistDeg2 = 1 / 2,
                tdistDeg4 = 1 / 4,
                rlang::abort("Invalid distribution type")
)

# Simulate s samples from distribution specified by type
set.seed(opt$seed)
samples <- purrr::map(1:s, \(i) rmvt(opt$n, sigma, 1 / gamma, mu))

# Create tibble of samples
d <- ncol(samples[[1]])
samples <- do.call(cbind, samples)
colnames(samples) <- paste0(c("x", "y"), rep(1:s, each = d))
samples <- tibble::as_tibble(samples)

# Calculate estimates
elliptical_estimates <- as.list(rep(NA, s))
depth_estimates <- as.list(rep(NA, s))
cli::cli_progress_bar("Compute estimates", total = s)
for (i in 1:s) {
  data <- samples %>%
    select(num_range(c("x", "y"), i)) %>%
    as.matrix
  est <- robustbase::covMcd(data, alpha = 0.5)
  
  e_i <- elliptical_extreme_qregion(data, est$cov, p, k, m)
  elliptical_estimates[[i]] <- e_i
  
  d_i <- depth_extreme_qregion(data, p, k, m)
  depth_estimates[[i]] <- d_i
  cli::cli_progress_update()
}

# Create tibble of estimates
column_names <- paste0(c("x", "y"), rep(1:s, each = d))

elliptical_estimates <- do.call(cbind, elliptical_estimates)
colnames(elliptical_estimates) <- column_names
elliptical_estimates <- tibble::as_tibble(elliptical_estimates)

depth_estimates <- do.call(cbind, depth_estimates)
colnames(depth_estimates) <- column_names
depth_estimates <- tibble::as_tibble(depth_estimates)

# Calculate theoretical quantile region
real <- tdist_extreme_region(sigma, gamma, p, m)
colnames(real) <- c("x", "y")
real <- tibble::as_tibble(real)

# Write data
filename <- paste0("type_", opt$type,
                   "_n_", opt$n,
                   "_seed_", opt$seed,
                   "_p_", opt$p,
                   "_k_", opt$k, ".csv"
                   )

readr::write_csv(elliptical_estimates,
                 paste0("data/elliptical-estimates/", filename)
                 )

readr::write_csv(depth_estimates,
                 paste0("data/depth-estimates/", filename)
                 )

filename <- paste0("type_", opt$type, "_p_", opt$p, ".csv")
readr::write_csv(real, paste0("data/real-regions/", filename))

filename <- paste0("type_", opt$type,
                   "_n_", opt$n,
                   "_seed_", opt$seed, ".csv"
                   )
readr::write_csv(samples, paste0("data/samples/", filename))
