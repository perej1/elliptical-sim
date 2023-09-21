# All the needed packages and helper functions
library(optparse)
library(mvtnorm)
library(ggplot2)
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(testthat))


#' Generate m equally spaced points from a circle
#'
#' @param m Integer, number of points.
#'
#' @return Double matrix of points, one row represents one point.
get_ball_mesh <- function(m) {
  w <- seq(0, 2 * pi, length.out = m)
  cbind(cos(w), sin(w))
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


#' Compute theoretical quantile region for t-distribution
#'
#' @param sigma Double matrix, Scatter matrix.
#' @param gamma Double, extreme value index.
#' @param p Double, probability mass in quantile region.
#' @param m Integer, number of points to return.
#'
#' @return Double matrix, m points from the boundary of the quantile region.
tdist_extreme_region <- function(sigma, gamma, p, m) {
  w <- get_ball_mesh(m)
  lambda <- sqrtmat(d * stats::qf(1 - p, d, 1 / gamma) * sigma)
  w %*% t(lambda)
}


#' Estimate elliptical extreme quantile region
#'
#' Consistent for heavy-tailed elliptical distributions under some other
#' technical assumptions.
#'
#' @param data Double matrix of observations, each row represents one
#'   observation.
#' @param mu_est Double vector, estimate of the location.
#' @param sigma_est Double matrix, estimate of the scatter.
#' @param p Double, probability in quantile region.
#' @param k Integer, threshold for the sample from the tail.
#' @param m Integer, number of points to return.
#'
#' @return Double matrix, m points from the boundary of the quantile region.
elliptical_extreme_qregion <- function(data, mu_est, sigma_est, p, k, m) {
  n <- nrow(data)
  w <- get_ball_mesh(m)

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
#' @param m Integer, number of points to return.
#'
#' @return Double matrix, m points from the boundary of the quantile region.
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
  nu_s_hat <- 1 / k *
    sum(radius / u_est >= hd_w_nu_hat_star[x_approx_w]^gamma_est)

  r <- u_est * (k * nu_s_hat / (n * p))^gamma_est * hd_w_nu_hat_star^gamma_est
  r * w
}


#' Calculate estimation error for elliptical extreme quantile region estimate
#'
#' Estimate \eqn{\mathbb{P}(X \in Q \triangle \hat Q)}, where \eqn{Q}
#' is the theoretical quantile region and \eqn{\hat Q} is the estimate.
#' Estimation of the integral \eqn{\mathbb{P}(X \in Q \triangle \hat Q)} is done
#' with Riemann sum in polar coordinates.
#'
#' @param real Double matrix, describes theoretical quantile region.
#' @param estimate Double matrix, describes estimated quantile region.
#' @param m1 Integer, discretization level for the angle.
#' @param m2 Integer, discretization level for the radius.
#' @param f density function.
#'
#' @return Double, estimation error.
compute_error <- function(real, estimate, m1, m2, f) {
  res <- 0

  ball <- get_ball_mesh(m1)
  theta <- rep(NA, m1)
  for (i in 1:m1) {
    theta[i] <- acos(ball[i, 1])
    if (ball[i, 2] < 0) {
      theta[i] <- 2 * pi - theta[i]
    }
  }

  r_real <- apply(real, 1, norm, type = "2")
  r_estimate <- apply(estimate, 1, norm, type = "2")

  ball_real <- sweep(real, 1, r_real, "/")
  ball_estimate <- sweep(estimate, 1, r_estimate, "/")

  ind_real <- apply(ball %*% t(ball_real), 1, which.max)
  ind_estimate <- apply(ball %*% t(ball_estimate), 1, which.max)

  r_real <- r_real[ind_real]
  r_estimate <- r_estimate[ind_estimate]

  for (i in 1:m1) {
    r_seq <- seq(min(r_estimate[i], r_real[i]),
                 max(r_estimate[i], r_real[i]),
                 length.out = m2)

    x <- r_seq * cos(theta[i])
    y <- r_seq * sin(theta[i])
    resi <- 0
    for (j in 1:m2) {
      resi <- resi + f(c(x[j], y[j])) * r_seq[j]
    }
    res <- res + abs(r_real[i] - r_estimate[i]) * resi
  }
  2 * pi / (m1 * m2) * res
}


#' Plot estimates and real region
#'
#' @param data Preprocessed data to plot, see summarise.R for preprocessing.
#'
#' @return The ggplot object.
plot_data <- function(data) {
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
