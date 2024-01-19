# Perform simulations for high dimensional scenarios. Output median errors
# corresponding certain combinations of d, n, p for optimal k and plot median
# errors for different dimensions.
source("functions.R")


#' Compute median conservative errors from s repetitions of a certain
#' simulation scenario
#'
#' @param s Integer, number of repetitions.
#' @param d Integer, number of dimensions.
#' @param n Integer, sample size.
#' @param p Character, p corresponding to quantile region, possible options are
#'   "low", "medium" and "high".
#' @param k Character, threshold for the sample from the tail, possible options
#'   are "small", "medium" and "large".
#' @param seed Integer, seed for set.seed().
#'
#' @return List of the median error and number of scenarios where the error
#' computation failed.
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

  k <- switch(k,
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


#' Save median errors for corresponding a certain parameter combination when
#'   dimensions run from 2 to 30
#'
#' @param n Integer, sample size.
#' @param k Character, threshold for the sample from the tail, possible options
#'   are "small", "medium" and "large".
#' @param p Character, p corresponding to quantile region, possible options are
#'   "low", "medium" and "high".
#'
#' @return Nothing, saves median errors and the corresponding plot.
plot_and_save <- function(n, k, p) {
  # Set parameter values
  seed <- 278
  s <- 100
  d <- 2:30

  # Computation of errors for each parameter combination
  errors <- expand.grid(s = s, d = d, n = n, p = p, k = k, seed = seed,
                        stringsAsFactors = FALSE) %>%
    rowwise() %>%
    mutate(median_err = list(simulate(s, d, n, p, k, seed))) %>%
    tidyr::unnest_wider(median_err)

  # Plotting
  g <- ggplot(errors, aes(d, error_med)) +
    geom_point() +
    geom_line() +
    ylim(0, NA) +
    xlab("Dimension") +
    ylab("Relative error")

  # Save data and the figure
  filename <- stringr::str_c("n_", n, "_k_", k, "_p_", p)
  readr::write_csv(errors,
                   stringr::str_c("high-dim-data/data/", filename, ".csv"))
  ggsave(stringr::str_c("high-dim-data/figures/", filename, ".jpg"),
         plot = g, width = 7, height = 7, dpi = 1000)
}

n <- c(1000, 5000)
k <- c("small", "medium", "large")
p <- c("low", "medium", "high")
args <- expand.grid(n = n, k = k, p = p, stringsAsFactors = FALSE)

cli::cli_progress_bar("High dimensional simulations ongoing...",
                      total = nrow(args))

for (i in seq_len(nrow(args))) {
  arg <- args[i, ]
  plot_and_save(arg$n, arg$k, arg$p)
  cli::cli_progress_update()
}
