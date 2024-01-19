# Unit tests for the function compute_error_elliptical
library(testthat)
source("functions.R")

alpha <- c(1, 5, 7)
p <- 0.05
d <- c(2, 5, 10)
p_wrong <- c(0.03, 0.07)

test_that("Error computed correctly when scatter is perfectly estimated", {
  for (alpha_i in alpha) {
    for (p_wrong_i in p_wrong) {
      for (d_i in d) {
        sigma_hat <- diag(d_i)
        r_hat <- sqrt(d_i * stats::qf(1 - p_wrong_i, d_i, alpha_i))
        err_est <- compute_error_elliptical(sigma_hat, r_hat, alpha_i, p)
        expect_equal(abs(p - p_wrong_i) / p, err_est)
      }
    }
  }
})

test_that("Scaling does not change results", {
  scale <- c(0.5, 2)
  for (alpha_i in alpha) {
    for (p_wrong_i in p_wrong) {
      for (d_i in d) {
        for (scale_i in scale) {
          sigma_hat <- diag(d_i)
          r_hat <- sqrt(d_i * stats::qf(1 - p_wrong_i, d_i, alpha_i))
          err_est <- compute_error_elliptical(scale_i^(-2) * sigma_hat,
                                              scale_i * r_hat, alpha_i, p)
          expect_equal(abs(p - p_wrong_i) / p, err_est)
        }
      }
    }
  }
})


test_that("Error behaves ok when scatter estimate is close to theoretical", {
  epsilon <- c(0.01, 0.05, 0.1)
  epsilon <- c(epsilon, -epsilon)

  for (alpha_i in alpha) {
    for (p_wrong_i in p_wrong) {
      for (d_i in d) {
        for (e_i in epsilon) {
          d_reverse <- diag(d_i)[d_i:1, ]
          sigma_hat <- diag(d_i) + e_i * d_reverse
          r_hat <- sqrt(d_i * stats::qf(1 - p_wrong_i, d_i, alpha_i))
          err_est <- compute_error_elliptical(sigma_hat, r_hat, alpha_i, p)
          expect_equal(abs(p - p_wrong_i) / p, err_est, tol = 0.4)
        }
      }
    }
  }
})
