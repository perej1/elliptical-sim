# Unit test function compute_error
source("functions.R")


set.seed(278)
d <- 2
m1 <- 1000
m2 <- 10
mu <- rep(0, 2)
p <- c(1 / 1000, 1 / 1500)
err <- p[1] - p[2]
gamma <- c(0.1, 0.5, 1)


test_that("Error computation works for spherical t-distriributions", {
  sigma <- diag(2)

  for (gamma_i in gamma) {
    real <- tdist_extreme_region(sigma, gamma_i, p[1], m1)
    estimate <- tdist_extreme_region(sigma, gamma_i, p[2], m1)
    err_est <- compute_error(real, estimate, m1, m2, gamma_i, sigma)
    expect_equal(err / p[1], err_est / p[1], tolerance = 0.01)
  }
})


test_that("Error computation works for elliptical t-distriributions", {
  sigma <- matrix(c(11, 10.5, 10.5, 11.25), byrow = TRUE, ncol = 2)

  for (gamma_i in gamma) {
    real <- tdist_extreme_region(sigma, gamma_i, p[1], m1)
    estimate <- tdist_extreme_region(sigma, gamma_i, p[2], m1)
    err_est <- compute_error(real, estimate, m1, m2, gamma_i, sigma)
    expect_equal(err / p[1], err_est / p[1], tolerance = 0.01)
  }
})
