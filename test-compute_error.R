# Unit test function compute_error
source("functions.R")


set.seed(278)
d <- 2
m1 <- 1000
m2 <- 100
mu <- rep(0, 2)
gamma <- c(0.25, 1)
tol <- 0.1


test_that("Error computation works for spherical t-distributions, when
          estimate is far away from real", {
  sigma <- diag(2)
  p <- c(1 / 1000, 1 / 100)
  err <- abs(p[1] - p[2])

  for (gamma_i in gamma) {
    real <- tdist_extreme_region(sigma, gamma_i, p[1], m1)
    estimate <- tdist_extreme_region(sigma, gamma_i, p[2], m1)
    err_est <- compute_error(real, estimate, m1, m2, gamma_i, sigma)
    print(paste0("Real: ", err / p[1]))
    print(paste0("Estimate: ", err_est / p[1]))
    expect_equal(err / p[1], err_est / p[1], tolerance = tol)
  }
})


test_that("Error computation works for elliptical t-distributions, when
          estimate is far away from real", {
  sigma <- matrix(c(11, 10.5, 10.5, 11.25), byrow = TRUE, ncol = 2)
  p <- c(1 / 1000, 1 / 100)
  err <- abs(p[1] - p[2])

  for (gamma_i in gamma) {
    real <- tdist_extreme_region(sigma, gamma_i, p[1], m1)
    estimate <- tdist_extreme_region(sigma, gamma_i, p[2], m1)
    err_est <- compute_error(real, estimate, m1, m2, gamma_i, sigma)
    print(paste0("Real: ", err / p[1]))
    print(paste0("Estimate: ", err_est / p[1]))
    expect_equal(err / p[1], err_est / p[1], tolerance = tol)
  }
})


test_that("Error computation works for spherical t-distributions when
          estimate is close to real", {
  sigma <- diag(2)
  p <- c(1 / 1000, 1 / 1500)
  err <- abs(p[1] - p[2])

  for (gamma_i in gamma) {
    real <- tdist_extreme_region(sigma, gamma_i, p[1], m1)
    estimate <- tdist_extreme_region(sigma, gamma_i, p[2], m1)
    err_est <- compute_error(real, estimate, m1, m2, gamma_i, sigma)
    print(paste0("Real: ", err / p[1]))
    print(paste0("Estimate: ", err_est / p[1]))
    expect_equal(err / p[1], err_est / p[1], tolerance = tol)
  }
})


test_that("Error computation works for elliptical t-distributions, when
          estimate is close to real", {
  sigma <- matrix(c(11, 10.5, 10.5, 11.25), byrow = TRUE, ncol = 2)
  p <- c(1 / 1000, 1 / 1500)
  err <- abs(p[1] - p[2])

  for (gamma_i in gamma) {
    real <- tdist_extreme_region(sigma, gamma_i, p[1], m1)
    estimate <- tdist_extreme_region(sigma, gamma_i, p[2], m1)
    err_est <- compute_error(real, estimate, m1, m2, gamma_i, sigma)
    print(paste0("Real: ", err / p[1]))
    print(paste0("Estimate: ", err_est / p[1]))
    expect_equal(err / p[1], err_est / p[1], tolerance = tol)
  }
})
