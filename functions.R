# All the needed packages and helper functions
library(ggplot2)
suppressPackageStartupMessages(library(dplyr))

#' Convert a spherical coordinate to a Cartesian coordinate
#'
#' @param radius Double, radius of the point.
#' @param theta Double vector, d - 1 angular coordinates.
#'
#' @return Double vector, d coordinates giving the location of the point in
#'   space.
spherical_to_cartesian <- function(radius, theta) {
  d <- length(theta) + 1
  cartesian <- rep(NA, d)
  cartesian[1] <- cos(theta[1])
  cartesian[d] <- prod(sin(theta))

  if (d > 2) {
    for (i in 2:(d - 1)) {
      cartesian[i] <- prod(sin(theta[1:(i - 1)])) * cos(theta[i])
    }
  }
  radius * cartesian
}


#' Generate m points from a (d - 1)-sphere
#'
#' @param d Integer, dimension of the ball, must be equal to or greater than 2.
#' @param m_angle Integer, number of points to return.
#'
#' @return List of two double matrices, the first gives the ball in spherical
#'   and the second in Cartesian coordinates. One row represents one point.
get_ball_mesh <- function(d, m_angle) {
  m <- ceiling(m_angle^(1 / (d - 1)))
  angle_list <- vector("list", d - 1)
  angle_list[[d - 1]] <- seq(0, 2 * pi, length.out = m)
  if (d > 2) {
    for (i in 1:(d - 2)) {
      angle_list[[i]] <- seq(0, pi, length.out = m)
    }
  }
  spherical <- as.matrix(expand.grid(angle_list))
  list(spherical = spherical,
       cartesian = t(apply(spherical, 1, spherical_to_cartesian, r = 1)),
       m_effective = nrow(spherical))
}


#' Compute square root of a positive definite matrix
#'
#' More precisely, function computes matrix lambda
#' s.t. sigma = lambda %*% lambda.
#'
#' @param sigma Double matrix, positive definite matrix.
#'
#' @return Double matrix, square root of the matrix.
sqrtmat <- function(sigma) {
  eigenval <- eigen(sigma)$values
  if (any(eigenval <= 0) || any(sigma != t(sigma))) {
    rlang::abort("`sigma` must be a symmetric positive definite matrix.")
  }
  eigenvec <- eigen(sigma)$vectors
  eigenvec %*% diag(eigenval^0.5) %*% t(eigenvec)
}


#' Compute theoretical (1-p)-quantile region for t-distribution
#'
#' @param sigma Double matrix, Scatter matrix.
#' @param gamma Double, extreme value index.
#' @param p Double, probability mass in quantile region.
#' @param m_angle Integer, number of points to return.
#'
#' @return Double matrix, m_angle points from the boundary of the quantile
#'   region.
tdist_extreme_region <- function(sigma, gamma, p, m_angle) {
  d <- ncol(sigma)
  w <- get_ball_mesh(d, m_angle)$cartesian
  lambda <- sqrtmat(d * stats::qf(1 - p, d, 1 / gamma) * sigma)
  w %*% t(lambda)
}


#' Estimate elliptical extreme quantile region
#'
#' Consistent for heavy-tailed elliptical distributions under some technical
#' assumptions.
#'
#' @param data Double matrix of observations, each row represents one
#'   observation.
#' @param mu_est Double vector, estimate of the location.
#' @param sigma_est Double matrix, estimate of the scatter.
#' @param p Double, probability in quantile region.
#' @param k Integer, threshold for the sample from the tail.
#' @param m_angle Integer, number of points to return.
#'
#' @return Double matrix, m_angle points from the boundary of the quantile
#'   region.
elliptical_extreme_qregion <- function(data, mu_est, sigma_est, p, k, m_angle) {
  n <- nrow(data)
  d <- ncol(data)
  w <- get_ball_mesh(d, m_angle)$cartesian

  # Center data
  data <- sweep(data, 2, mu_est, "-")

  # Approximate generating variate
  radius <- sqrt(stats::mahalanobis(data, FALSE, sigma_est, inverted = FALSE))
  radius_sort <- sort(radius, decreasing = FALSE)

  # Estimate extreme value index
  gamma_est <- mean((log(radius_sort[(n - k):n]) - log(radius_sort[n - k]))[-1])

  # Estimate extreme quantile of generating variate
  r_hat <- radius_sort[n - k] * (k / (n * p))^gamma_est

  # Estimate extreme quantile region
  lambda <- sqrtmat(r_hat^2 * sigma_est)
  sweep(w %*% t(lambda), 2, mu_est, "+")
}


#' Estimate halfspace depth extreme quantile region
#'
#' Consistent for multivariate regularly varying distributions under some other
#' technical assumptions. See "Y. He, J. H. Einmahl, Estimation of extreme
#' depth-based quantile regions" for more details.
#'
#' @param data Double matrix of observations, each row represents one
#'   observation.
#' @param p Double, probability in quantile region.
#' @param k Integer, threshold for the sample from the tail.
#' @param m_angle Integer, number of points to return.
#'
#' @return Double matrix, m_angle points from the boundary of the quantile
#'   region.
depth_extreme_qregion <- function(data, p, k, m_angle) {
  n <- nrow(data)
  d <- ncol(data)
  w <- get_ball_mesh(d, m_angle)$cartesian

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
  nu_s_hat <- 1 / k *
    sum(radius / u_est >= hd_w_nu_hat_star[x_approx_w]^gamma_est)

  r <- u_est * (k * nu_s_hat / (n * p))^gamma_est * hd_w_nu_hat_star^gamma_est
  r * w
}


#' Compute density contour at level beta for skewed t-distribution
#'
#' Computes points from the boundary of \eqn{\{x : f(x) \leq \beta\}}.
#'
#' @param mu Double vector, location.
#' @param sigma Double matrix, scatter matrix.
#' @param gamma Double, extreme value index.
#' @param alpha Double vector, Skewness parameter for the skewed t-distribution.
#' @param beta Double, level corresponding to density contour.
#' @param m_angle Integer, number of points to return.
#'
#' @return Double matrix, m_angle points from the density contour.
skew_t_contour_beta <- function(mu, sigma, gamma, alpha, beta, m_angle) {
  theta <- seq(0, 2 * pi, length.out = m_angle)
  r <- rep(NA, m_angle)
  for (i in 1:m_angle) {
    g <- function(rr) {
      x <- rr * cos(theta[i])
      y <- rr * sin(theta[i])
      sn::dmst(c(x, y), mu, sigma, alpha, 1 / gamma, log = FALSE) - beta
    }
    r[i] <- pracma::fzero(g, c(0, 1000))$x
  }

  x <- r * cos(theta)
  y <- r * sin(theta)
  cbind(x, y)
}


#' Compute theoretical (1-p)-quantile region for skewed t-distribution
#'
#' @param mu Double vector, location.
#' @param sigma Double matrix, scatter matrix.
#' @param gamma Double, extreme value index.
#' @param alpha Double vector, Skewness parameter for the skewed t-distribution.
#' @param p Double, probability mass in quantile region.
#' @param m_angle Integer, number of points to return.
#'
#' @return Double matrix, m_angle points from the boundary of the quantile
#'   region.
skew_t_contour_p <- function(mu, sigma, gamma, alpha, p, m_angle) {
  n <- 2 * ceiling(1 / p)
  f_sample <- apply(sn::rmst(n, mu, sigma, alpha, 1 / gamma), 1,
                    function(x) sn::dmst(x, mu, sigma, alpha, 1 / gamma, log = FALSE))
  beta <- quantile(f_sample, p)
  skew_t_contour_beta(mu, sigma, gamma, alpha, beta, m_angle)
}


#' Compute estimation error for elliptical extreme quantile region estimate
#'
#' Estimate \eqn{\mathbb{P}(X \in Q \triangle \hat Q)}, where \eqn{Q}
#' is the theoretical quantile region and \eqn{\hat Q} is the estimate.
#' Estimation of the integral \eqn{\mathbb{P}(X \in Q \triangle \hat Q)} is done
#' with Riemann sum in polar coordinates.
#'
#' @param real Double matrix, describes theoretical quantile region.
#' @param estimate Double matrix, describes estimated quantile region.
#' @param m_radius Integer, discretization level for the radius.
#' @param f Density function.
#'
#' @return Double, estimation error.
compute_error <- function(real, estimate, m_radius, f) {
  m_angle <- nrow(real)
  d <- ncol(real)
  ball <- get_ball_mesh(d, m_angle)
  if (m_angle != nrow(estimate)) {
    rlang::abort("`real` and `estimate` must have the number of rows")
  }
  res <- 0

  r_real <- apply(real, 1, norm, type = "2")
  r_estimate <- apply(estimate, 1, norm, type = "2")

  ball_real <- sweep(real, 1, r_real, "/")
  ball_estimate <- sweep(estimate, 1, r_estimate, "/")

  for (i in 1:m_angle) {
    theta <- ball$spherical[i, ]
    i_real <- which.max(ball_real %*% ball$cartesian[i, ])
    i_estimate <- which.max(ball_estimate %*% ball$cartesian[i, ])
    r_seq <- seq(min(r_estimate[i_estimate], r_real[i_real]),
                 max(r_estimate[i_estimate], r_real[i_real]),
                 length.out = m_radius)

    cartesian <- t(sapply(r_seq, spherical_to_cartesian, theta = theta))
    jdet_angle <- ifelse(d > 2, prod(sin(theta[1:(d - 2)])^((d - 2):1)), 1)
    jdet <- r_seq^(d - 1) * jdet_angle
    resi <- sum(apply(cartesian, 1, f) * jdet)
    res <- res + abs(r_real[i_real] - r_estimate[i_estimate]) * resi
  }
  2 * pi^(d - 1) / (m_angle * m_radius) * sum(res)
}


#' Plot estimates and real region together for one setting
#'
#' @param real m_angle by 2 tibble, theoretical quantile region.
#' @param ellipse m_angle by 2 tibble, elliptical estimate.
#' @param depth m_angle by 2 tibble, depth based estimate.
#'
#' @return The ggplot object.
plot_real_estimate <- function(real, ellipse, depth) {
  m_angle <- nrow(real)
  data <- bind_rows(real, ellipse, depth) %>%
    mutate(group = rep(c("real", "ellipse", "depth"), each = m_angle))
  
  ggplot(data, aes(x = x, y = y)) +
    geom_path(aes(group = group, linetype = group), show.legend = FALSE) +
    coord_fixed() +
    theme(axis.title = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          legend.background = element_blank(),
          legend.key = element_blank(),
          legend.key.size = unit(1, "cm"),
          legend.text = element_text(size = 15),
          axis.text = element_text(size = 15)) +
    scale_linetype_manual(name = NULL,
                          values = c("real" = "solid", "ellipse" = "dashed",
                                     "depth" = "dotted"))
}


#' Plot estimate and theoretical quantile regions for different values of p.
#'
#' @param sample tibble of observations, one row represent one observation.
#' @param real_list List of tibbles, one element of the list represents
#'   theoretical quantile region corresponding a specific probability.
#' @param estimate_list List of tibbles, one element of the list represents
#'   estimated quantile region corresponding a specific probability.
#' @param p Vector of strings, labels for different values of probabilities
#'   corresponding to quantile regions.
#'
#' @return The ggplot object.
plot_sample_and_estimates <- function(sample, real_list, estimate_list, p) {
  m_angle <- nrow(real_list[[1]])
  estimates <- bind_rows(estimate_list) %>%
    mutate(group = rep(paste0("estimate_", p), each = m_angle))
  
  reals <- bind_rows(real_list) %>%
    mutate(group = rep(paste0("real_", p), each = m_angle))
  
  linetypes <- rep(c("solid", "dashed"), each = length(p))
  names(linetypes) <- c(paste0("real_", p), paste0("estimate_", p))
  
  g <- ggplot(as_tibble(sample), aes(x = x, y = y)) +
    geom_point() +
    geom_path(data = bind_rows(estimates, reals),
              aes(x = x, y = y, group = group, linetype = group),
              show.legend = FALSE) +
    coord_fixed() +
    theme(axis.title = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          legend.background = element_blank(),
          legend.key = element_blank(),
          legend.key.size = unit(1, "cm"),
          legend.text = element_text(size = 15),
          axis.text = element_text(size = 15)) +
    scale_linetype_manual(name = NULL,
                          values = linetypes)
}


#' Boxplot errors for elliptical and depth based estimator for one scenario.
#'
#' @param elliptical s times one tibble, errors for elliptical estimator for s
#'   repetitions of one simulation scenario.
#' @param depth s times one tibble, errors for depth based estimator for s
#'   repetitions of one simulation scenario.
#'
#' @return The ggplot object.
boxplot_errors <- function(elliptical, depth) {
  s <- nrow(elliptical)
  elliptical <- rename(elliptical, err = 1)
  depth <- rename(depth, err = 1)
  data <- rbind(elliptical, depth) %>%
    mutate(group = ordered(rep(c("Elliptical", "Depth"), each = s),
                           levels = c("Elliptical", "Depth")))
  g <- ggplot(data, aes(x = group, y = err)) +
    geom_boxplot() +
    theme(axis.title = element_blank(),
          panel.background = element_blank(),
          legend.background = element_blank(),
          legend.key = element_blank(),
          legend.key.size = unit(1, "cm"),
          axis.text = element_text(size = 15))
}
