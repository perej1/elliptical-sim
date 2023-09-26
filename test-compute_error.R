# Unit tests for the function compute_error
source("functions.R")


m_angle <- c(100, 10000)
m_radius <- c(100, 100)
mu <- list(rep(0, 2), rep(0, 3))
sigma <- list(
  matrix(c(11, 10.5, 10.5, 11.25), byrow = TRUE, ncol = 2),
  matrix(c(8, 7.5, -2.25, 7.5, 15, 0.45, -2.25, 0.45, 2),
         byrow = TRUE, ncol = 3))
gamma <- c(0.25, 1)
gamma_i <- 0.25
tol <- 0.1


test_that("Volume of a sphere is calculated correctly", {
  for (d in 2:3) {
    ball1 <- get_ball_mesh(d, m_angle[d - 1])$cartesian
    ball2 <- 1.1 * ball1
    err <- (1.1^d - 1) * pi^(d / 2) / gamma(d / 2 + 1)
    tic()
    err_est <- compute_error(ball1, ball2, m_radius[d - 1], function(x) 1)
    toc()
    print(paste0("dim = ", d))
    print(paste0("Real: ", err))
    print(paste0("Estimate: ", err_est))
    expect_equal(err, err_est, tolerance = tol)
  }
})


test_that("Error computation works for spherical t-distributions, when
          estimate is far away from real", {
  for (d in 2:3) {
    p <- c(1 / 1000, 1 / 100)
    err <- abs(p[1] - p[2])
    
    for (gamma_i in gamma) {
      f <- function(x) dmvt(x, mu[[d - 1]], diag(d), 1 / gamma_i, log = FALSE)
      real <- tdist_extreme_region(diag(d), gamma_i, p[1], m_angle[d - 1])
      estimate <- tdist_extreme_region(diag(d), gamma_i, p[2], m_angle[d - 1])
      tic()
      err_est <- compute_error(real, estimate, m_radius[d - 1], f)
      toc()
      print(paste0("dim = ", d, ", gamma = ", gamma_i))
      print(paste0("Real: ", err / p[1]))
      print(paste0("Estimate: ", err_est / p[1]))
      expect_equal(err / p[1], err_est / p[1], tolerance = tol)
    }
  }
})


test_that("Error computation works for elliptical t-distributions, when
          estimate is far away from real", {
  for (d in 2:3) {
    p <- c(1 / 1000, 1 / 100)
    err <- abs(p[1] - p[2])
    
    for (gamma_i in gamma) {
      f <- function(x) dmvt(x, mu[[d - 1]], sigma[[d - 1]], 1 / gamma_i, log = FALSE)
      real <- tdist_extreme_region(sigma[[d - 1]], gamma_i, p[1], m_angle[d - 1])
      estimate <- tdist_extreme_region(sigma[[d - 1]], gamma_i, p[2], m_angle[d - 1])
      tic()
      err_est <- compute_error(real, estimate, m_radius[d - 1], f)
      toc()
      print(paste0("dim = ", d, ", gamma = ", gamma_i))
      print(paste0("Real: ", err / p[1]))
      print(paste0("Estimate: ", err_est / p[1]))
      expect_equal(err / p[1], err_est / p[1], tolerance = tol)
    }
  }
})


test_that("Error computation works for spherical t-distributions, when
          estimate is close to real", {
  for (d in 2:3) {
    p <- c(1 / 1000, 1 / 1500)
    err <- abs(p[1] - p[2])
    
    for (gamma_i in gamma) {
      f <- function(x) dmvt(x, mu[[d - 1]], diag(d), 1 / gamma_i, log = FALSE)
      real <- tdist_extreme_region(diag(d), gamma_i, p[1], m_angle[d - 1])
      estimate <- tdist_extreme_region(diag(d), gamma_i, p[2], m_angle[d - 1])
      tic()
      err_est <- compute_error(real, estimate, m_radius[d - 1], f)
      toc()
      print(paste0("dim = ", d, ", gamma = ", gamma_i))
      print(paste0("Real: ", err / p[1]))
      print(paste0("Estimate: ", err_est / p[1]))
      expect_equal(err / p[1], err_est / p[1], tolerance = tol)
    }
  }
})


test_that("Error computation works for elliptical t-distributions, when
          estimate is close to real", {
  for (d in 2:3) {
    p <- c(1 / 1000, 1 / 1500)
    err <- abs(p[1] - p[2])
    
    for (gamma_i in gamma) {
      f <- function(x) dmvt(x, mu[[d - 1]], sigma[[d - 1]], 1 / gamma_i, log = FALSE)
      real <- tdist_extreme_region(sigma[[d - 1]], gamma_i, p[1], m_angle[d - 1])
      estimate <- tdist_extreme_region(sigma[[d - 1]], gamma_i, p[2], m_angle[d - 1])
      tic()
      err_est <- compute_error(real, estimate, m_radius[d - 1], f)
      toc()
      print(paste0("dim = ", d, ", gamma = ", gamma_i))
      print(paste0("Real: ", err / p[1]))
      print(paste0("Estimate: ", err_est / p[1]))
      expect_equal(err / p[1], err_est / p[1], tolerance = tol)
    }
  }
})
