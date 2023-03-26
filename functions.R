library(optparse)
library(mvtnorm)
library(ggplot2)
suppressPackageStartupMessages(library(dplyr))


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
#' @param p Double, probability in quantile region.
#' @param m Integer, number of points to return.
#'
#' @return Double matrix, m points from the boundary of quantile region.
tdist_extreme_region <- function(sigma, gamma, p, m) {
  w <- get_ball_mesh(m)
  lambda <- sqrtmat(d * stats::qf(1 - p, d, 1 / gamma) * sigma)
  w %*% t(lambda)
}


#' Estimate elliptical extreme quantile region.
#' 
#' Consistent for heavy-tailed distributions under some other technical
#' assumptions.
#'
#' @param data Double matrix of observations, each row represents one
#'   observation.
#' @param mu_est Double vector, estimate of location.
#' @param sigma_est Double matrix, estimate of scatter.
#' @param p Double, probability in quantile region.
#' @param k Integer, threshold for the sample from the tail.
#' @param m Integer, number of points to return.
#'
#' @return Double matrix, m points from the boundary of quantile region.
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


density_clover <- function(x, y) {
  r_0 <- 1.2481
  if (x^2 + y^2 < r_0) {
    a <- 3 * r_0^4 * (1 + r_0^6)^(-3/2) / (10 * pi)
    up <- 4 * (x^2 + y^2)^2 - 32 * x^2 * y^2
    down <- r_0 * (x^2 + y^2)^(3/2)
    a * (5 + up / down)
  } else {
    up <- 3 * (9 * (x^2 + y^2)^2 - 32 * x^2 * y^2)
    down <- 10 * pi * (1 + (x^2 + y^2)^3)^(3/2)
    up / down
  }
}


gen_clover <- function(n) {
  xy_max <- 30
  z_max <- density_clover(0.0001, 0.0001)
  ret <- matrix(NA, nrow = n, ncol = 2)
  for (i in 1:n) {
    x <- runif(1, min = -xy_max, max = xy_max)
    y <- runif(1, min = -xy_max, max = xy_max)
    z <- runif(1, min = 0, max = z_max)
    while (z > density_clover(x, y)) {
      x <- runif(1, min = -xy_max, max = xy_max)
      y <- runif(1, min = -xy_max, max = xy_max)
      z <- runif(1, min = 0, max = z_max)
    }
    ret[i, ] <- c(x,y)
  }
  ret
}


clover_contour_beta <- function(beta, m) {
  theta <- seq(0, 2*pi, length.out = m)
  r <- rep(NA, m)
  for (i in 1:m) {
    g <- function(x) {
      r_0 <- 1.2481
      if (x < r_0) {
        a <- 3 * r_0^4 * (1 + r_0^6)^(-3/2) / (10 * pi)
        5 * a * r_0 + 4 * x * a * (1 - 2 * (sin(2 * theta[i]))^2) - beta * r_0
      } else {
        a <- 3 * (9 - 8 * (sin(2 * theta[i]))^2) / (10 * pi)
        a * x^4 - beta * (1 + x^6)^(3/2)
      }
    }
    r[i] <- pracma::fzero(g, c(0, 1000))$x
  }
  
  x <- r * cos(theta)
  y <- r * sin(theta)
  cbind(x, y)
}


clover_contour_p <- function(p, m) {
  n <- 2 * ceiling(1 / p)
  f_sample <- apply(gen_clover(n), 1, function (x) density_clover(x[1], x[2]))
  beta <- quantile(f_sample, p)
  clover_contour_beta(beta, m)
}


compute_error <- function(real, estimate, m1, m2, f, sigma) {
  l <- solve(sqrtmat(sigma))
  estimate <- estimate %*% t(l)
  real <- real %*% t(l)
  res <- 0
  coord <- matrix(NA, ncol = 2, nrow = m1 * m2)
  
  r_real <- apply(real, 1, norm, type = "2")
  r_estimate <- apply(estimate, 1, norm, type = "2")
  
  ball_real <- sweep(real, 1, r_real, "/")
  ball_estimate <- sweep(estimate, 1, r_estimate, "/")
  
  theta <- rep(NA, m1)
  for (i in 1:m1) {
    theta[i] <- acos(ball_real[i, 1])
    if (ball_real[i, 2] < 0) {
      theta[i] <- 2 * pi - theta[i]
    }
  }
  
  ind <- apply(ball_real %*% t(ball_estimate), 1, which.max)
  r_estimate <- r_estimate[ind]
  
  for (i in 1:m1) {
    r_seq <- seq(min(r_estimate[i], r_real[i]),
                 max(r_estimate[i], r_real[i]),
                 length.out = m2)
    
    x <- r_seq * cos(theta[i])
    y <- r_seq * sin(theta[i])
    coord[((i-1)*m2+1):(i*m2), ] <- cbind(x,y)
    resi <- 0
    for (j in 1:m2) {
      resi <- resi + f(c(x[j], y[j])) * r_seq[j]
    }
    res <- res + abs(r_real[i] - r_estimate[i]) * resi
  }
  list(res = 2 * pi / (m1 * m2) * res, coord = coord)
}
