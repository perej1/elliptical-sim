# Perform simulation for one scenario.
source("functions.R")

option_list <- list(
  make_option("--type", type = "character", default = "cauchyAff",
              help = "Distribution type"),
  make_option("--n", type = "integer", default = 500,
              help = "Sample size"),
  make_option("--p", type = "character", default = "low",
              help = "Probability mass outside quantile region"),
  make_option("--k", type = "character", default = "small",
              help = "Sample size of the tail"),
  make_option("--seed", type = "integer", default = 278,
              help = "Set seed for sampling")
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
  
  e_i <- elliptical_extreme_qregion(data, mu, est$cov, p, k, m)
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

# Density function
f <- function (x) {
  dmvt(x, mu, diag(2), 1 / gamma, log = FALSE)
}

# Calculate errors
elliptical_err <- rep(NA, s)
depth_err <- rep(NA, s)
cli::cli_progress_bar("Compute errors", total = s)
for (i in 1:s) {
  e_est <- elliptical_estimates %>%
    select(num_range(c("x", "y"), i)) %>%
    as.matrix
  d_est <- depth_estimates %>%
    select(num_range(c("x", "y"), i)) %>%
    as.matrix
  elliptical_err[i] <- compute_error(as.matrix(real), e_est, m, 10, f, sigma)$res / p
  depth_err[i] <- tryCatch(compute_error(as.matrix(real), d_est, m, 10, f, sigma)$res,
                           error = function(err) NA) / p
  cli::cli_progress_update()
}

errors <- tibble(elliptical = elliptical_err, depth = depth_err)

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

readr::write_csv(errors,
                 paste0("data/errors/", filename)
)

filename <- paste0("type_", opt$type, "_p_", opt$p, ".csv")
readr::write_csv(real, paste0("data/real-regions/", filename))

filename <- paste0("type_", opt$type,
                   "_n_", opt$n,
                   "_seed_", opt$seed, ".csv"
                   )
readr::write_csv(samples, paste0("data/samples/", filename))
