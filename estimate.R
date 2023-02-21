library(mvtnorm)
library(tibble)
library(rlang)
library(robustbase)

get_ball_mesh <- function(d, m) {
  if (d == 2) {
    w <- seq(0, 2 * pi, length.out = m)
    cbind(cos(w), sin(w))
  } else if (d == 3) {
    inc <- seq(0, pi, length.out = m)
    azi <- seq(0, 2 * pi, length.out = m)
    cbind(sin(inc) * cos(azi), sin(inc) * sin(azi), cos(inc))
  } else {
    abort("Dimensions of data should be equal to two or three.")
  }
}

sqrtmat <- function(sigma) {
  eigenval <- eigen(sigma)$values
  if (any(eigenval <= 0) || any(sigma != t(sigma))) {
    abort("`sigma` must be a symmetric positive definite matrix.")
  }
  eigenvec <- eigen(sigma)$vectors
  eigenvec %*% diag(eigenval^0.5) %*% t(eigenvec)
}

depth_extreme_qregion <- function(data, p, k, m) {
  d <- ncol(data)
  n <- nrow(data)
  w <- get_ball_mesh(d, m)
  
  # Center the data
  center_est <- covMcd(data, alpha = 0.5)$center
  data <- sweep(data, 2, center_est, "-")
  
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
  coord <- sweep(r * w, 2, center_est, "+")
  if (d == 2) {
    tibble(x = coord[, 1], y = coord[, 2])
  }
  else if (d == 3) {
    tibble(x = coord[, 1], y = coord[, 2], z = coord[, 3])
  }
}

elliptical_extreme_qregion <- function(data, p, k, m) {
  d <- ncol(data)
  n <- nrow(data)
  w <- get_ball_mesh(d, m)
  est <- covMcd(data, alpha = 0.5)
  
  # Center data
  data <- sweep(data, 2, est$center, "-")
  
  # Approximate generating variate
  radius <- sqrt(mahalanobis(data, FALSE, est$cov, inverted = FALSE))
  radius_sort <- sort(radius, decreasing = FALSE)
  
  gamma_est <- mean((log(radius_sort[(n - k):n]) - log(radius_sort[n - k]))[-1])
  lambda <- sqrtmat((radius_sort[n - k] * (k / (n * p))^gamma_est)^2 * est$cov)
  
  coord <- sweep(w %*% t(lambda), 2, est$center, "+")
  if (d == 2) {
    tibble(x = coord[, 1], y = coord[, 2])
  }
  else if (d == 3) {
    tibble(x = coord[, 1], y = coord[, 2], z = coord[, 3])
  }
}

set.seed(123)
n <- 5000
k <- 400
m <- 1000
p <- 1 / 5000
mu <- c(10000, 10000)
sigma <- matrix(c(1, 0.9, 0.9, 2), byrow = TRUE, ncol = 2)
data <- rmvt(n, sigma, df = 1, delta = mu)

ret_depth <- depth_extreme_qregion(data, p, k, m)
ret_elliptical <- elliptical_extreme_qregion(data, p, k, m)

plot(ret_depth$x, ret_depth$y, type = "l")
points(ret_elliptical$x, ret_elliptical$y, type = "l")
points(data[, 1], data[, 2])
