library(optparse)
suppressPackageStartupMessages(library(dplyr))

get_ball_mesh <- function(d, m) {
  if (d == 2) {
    w <- seq(0, 2 * pi, length.out = m)
    cbind(cos(w), sin(w))
  } else if (d == 3) {
    inc <- seq(0, pi, length.out = m)
    azi <- seq(0, 2 * pi, length.out = m)
    cbind(sin(inc) * cos(azi), sin(inc) * sin(azi), cos(inc))
  } else {
    rlang::abort("Dimensions of data should be equal to two or three.")
  }
}

sqrtmat <- function(sigma) {
  eigenval <- eigen(sigma)$values
  if (any(eigenval <= 0) || any(sigma != t(sigma))) {
    rlang::abort("`sigma` must be a symmetric positive definite matrix.")
  }
  eigenvec <- eigen(sigma)$vectors
  eigenvec %*% diag(eigenval^0.5) %*% t(eigenvec)
}

tdist_extreme_region <- function(mu, sigma, gamma, p, m) {
  d <- length(mu)
  w <- get_ball_mesh(d, m)
  lambda <- sqrtmat(d * stats::qf(1 - p, d, 1 / gamma) * sigma)
  sweep(w %*% t(lambda), 2, mu, "+")
}

elliptical_extreme_qregion <- function(data, mu_est, sigma_est, p, k, m) {
  d <- ncol(data)
  n <- nrow(data)
  w <- get_ball_mesh(d, m)
  
  # Center data
  data <- sweep(data, 2, mu_est, "-")
  
  # Approximate generating variate
  radius <- sqrt(stats::mahalanobis(data, FALSE, sigma_est, inverted = FALSE))
  radius_sort <- sort(radius, decreasing = FALSE)
  
  gamma_est <- mean((log(radius_sort[(n - k):n]) - log(radius_sort[n - k]))[-1])
  lambda <- sqrtmat((radius_sort[n - k] * (k / (n * p))^gamma_est)^2 * sigma_est)
  
  sweep(w %*% t(lambda), 2, mu_est, "+")
}

depth_extreme_qregion <- function(data, mu_est, p, k, m) {
  d <- ncol(data)
  n <- nrow(data)
  w <- get_ball_mesh(d, m)
  
  # Center the data
  data <- sweep(data, 2, mu_est, "-")
  
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
  sweep(r * w, 2, mu_est, "+")
}

# Argument list
option_list <- list(
  make_option("--type", type = "character", default = "cauchy",
              help = "Distribution type"),
  make_option("--n", type = "integer", default = 100,
              help = "Sample size"),
  make_option("--p", type = "character", default = "low",
              help = "Probability mass outside quantile region"),
  make_option("--k", type = "character", default = "large",
              help = "Sample size of the tail")
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Global constants
d <- ifelse(opt$type == "cauchy3d", 3, 2)
m <- 1000
s <- 100

# Read data
filename <- paste0(opt$type, "_", opt$n, ".csv")
samples <- readr::read_csv(paste0("data/samples/", filename), show_col_types = FALSE)

# Set values for p, k, m and extreme value index
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

if (opt$type == "cauchyAff") {
  mu <- c(100, -250)
  sigma <- matrix(c(11, 10.5, 10.5, 11.25), byrow = TRUE, ncol = 2)
} else {
  mu <- rep(0, d)
  sigma <- diag(d)
}

# Calculate estimates
elliptical_estimates <- as.list(rep(NA, s))
depth_estimates <- as.list(rep(NA, s))

for (i in 1:s) {
  data <- samples %>%
    select(num_range(c("x", "y", "z")[1:d], i)) %>%
    as.matrix
  est <- robustbase::covMcd(data, alpha = 0.5)
  
  e_i <- elliptical_extreme_qregion(data, est$center, est$cov, p, k, m)
  elliptical_estimates[[i]] <- e_i
  
  d_i <- depth_extreme_qregion(data, est$center, p, k, m)
  depth_estimates[[i]] <- d_i
}

# Create tibble of estimates
column_names <- paste0(c("x", "y", "z")[1:d], rep(1:s, each = d))

elliptical_estimates <- do.call(cbind, elliptical_estimates)
colnames(elliptical_estimates) <- column_names
elliptical_estimates <- tibble::as_tibble(elliptical_estimates)

depth_estimates <- do.call(cbind, depth_estimates)
colnames(depth_estimates) <- column_names
depth_estimates <- tibble::as_tibble(depth_estimates)

# Calculate theoretical quantile region
real <- tdist_extreme_region(mu, sigma, gamma, p, m)
colnames(real) <- c("x", "y", "z")[1:d]
real <- tibble::as_tibble(real)

# Write data
filename <- paste0(opt$type, "_", opt$n, "_", opt$p, "_", opt$k, ".csv")
readr::write_csv(elliptical_estimates, paste0("data/elliptical-estimates/", filename))
readr::write_csv(depth_estimates, paste0("data/depth-estimates/", filename))

filename <- paste0(opt$type, "_", opt$n, "_", opt$p, ".csv")
readr::write_csv(real, paste0("data/real-regions/", filename))
