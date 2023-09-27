# Perform simulation for one scenario.
source("functions.R")
tictoc::tic()

option_list <- list(
  make_option("--type", type = "character", default = "tdistDeg4",
              help = "Distribution type"),
  make_option("--s", type = "integer", default = 1,
              help = "Dimensions"),
  make_option("--d", type = "integer", default = 2,
              help = "Dimensions"),
  make_option("--n", type = "integer", default = 1000,
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
m_radius <- 100
m_angle <- c(100, 10000)
sigma_list <- list(
  matrix(c(11, 10.5, 10.5, 11.25), byrow = TRUE, ncol = 2),
  matrix(c(8, 7.5, -2.25, 7.5, 15, 0.45, -2.25, 0.45, 2),
         byrow = TRUE, ncol = 3)
)
alpha <- c(1, -3)

# Scatter and location
if (opt$type == "cauchyAff") {
  sigma <- sigma_list[[opt$d - 1]]
} else {
  sigma <- diag(opt$d)
}
mu <- rep(0, opt$d)

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
  tdistDeg4 = 1 / 4,
  tdistSkew = 1 / 4,
  rlang::abort("Invalid distribution type")
)

# Compute theoretical quantile region and set density function
if (opt$type == "tdistSkew") {
  f <- function(x) sn::dmst(x, mu, sigma, alpha, nu = 1 / gamma)
  real <- NA
} else {
  f <- function(x) mvtnorm::dmvt(x, mu, sigma, df = 1 / gamma, log = FALSE)
  real <- tdist_extreme_region(sigma, gamma, p, m_angle[opt$d - 1])
}

# One round of simulations
simulate <- function(i) {
  # Simulate sample
  if (opt$type == "tdistSkew") {
    data <- sn::rmst(opt$n, mu, sigma, alpha, 1 / gamma)
  } else {
    data <- mvtnorm::rmvt(opt$n, sigma, 1 / gamma, mu)
  }

  # Estimation
  est <- robustbase::covMcd(data, alpha = 0.5)
  e_i <- elliptical_extreme_qregion(data, mu, est$cov, p, k, m_angle[opt$d - 1])
  d_i <- depth_extreme_qregion(data, p, k, m_angle[opt$d - 1])

  # Compute error
  e_err_i <- compute_error(real, e_i, m_radius, f) / p
  d_err_i <- tryCatch(compute_error(real, d_i, m_radius, f),
                      error = function(err) NA) / p

  list(samples = data,
    elliptical_estimates = e_i,
    depth_estimates = d_i,
    elliptical_err = e_err_i,
    depth_err = d_err_i
  )
}


set.seed(opt$seed)
res <- purrr::map(1:opt$s, ~simulate(.x)) %>%
  purrr::transpose()

coord <- c("x", "y", "z")[1:opt$d]

# Collect all samples in the same tibble
samples <- do.call(cbind, res$samples)
colnames(samples) <- paste0(coord, rep(1:opt$s, each = opt$d))
samples <- tibble::as_tibble(samples)

# Collect all the estimates in the same tibble
column_names <- paste0(coord, rep(1:opt$s, each = opt$d))

elliptical_estimates <- do.call(cbind, res$elliptical_estimates)
colnames(elliptical_estimates) <- column_names
elliptical_estimates <- tibble::as_tibble(elliptical_estimates)

depth_estimates <- do.call(cbind, res$depth_estimates)
colnames(depth_estimates) <- column_names
depth_estimates <- tibble::as_tibble(depth_estimates)

# Collect all errors in the same tibble
errors <- tibble(elliptical = purrr::flatten_dbl(res$elliptical_err),
                 depth = purrr::flatten_dbl(res$depth_err))

# Theoretical quantile region to tibble
colnames(real) <- coord
real <- tibble::as_tibble(real)

# Write data
filename <- paste0("type_", opt$type,
                   "_s_", opt$s,
                   "_d_", opt$d,
                   "_n_", opt$n,
                   "_p_", opt$p,
                   "_k_", opt$k,
                   "_seed_", opt$seed, ".csv")

readr::write_csv(elliptical_estimates,
                 paste0("sim-data/elliptical-estimates/", filename))

readr::write_csv(depth_estimates,
                 paste0("sim-data/depth-estimates/", filename))

readr::write_csv(errors,
                 paste0("sim-data/errors/", filename))

filename <- paste0("type_", opt$type,
                   "_d_", opt$d,
                   "_n_", opt$n,
                   "_p_", opt$p, ".csv")
readr::write_csv(real, paste0("sim-data/real-regions/", filename))

filename <- paste0("type_", opt$type,
                   "_s_", opt$s,
                   "_d_", opt$d,
                   "_n_", opt$n,
                   "_seed_", opt$seed, ".csv")
readr::write_csv(samples, paste0("sim-data/samples/", filename))
tictoc::toc()
