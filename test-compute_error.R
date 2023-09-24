# Unit tests for the function compute_error
source("functions.R")


set.seed(278)
m_angle <- 1000
m_radius <- 3
mu <- rep(0, 2)
gamma <- c(0.25, 1)
tol <- 0.1


test_that("Volume of a sphere is calculated correctly", {
  for (d in 2:3) {
    ball1 <- get_ball_mesh(d, m_angle)$cartesian
    ball2 <- 1.1 * ball1
    err <- (1.1^d - 1) * pi^(d / 2) / gamma(d / 2 + 1)
    err_est <- compute_error(ball1, ball2, m_angle, m_radius, function(x) 1)
    print(paste0("Real: ", err))
    print(paste0("Estimate: ", err_est))
    expect_equal(err, err_est, tolerance = tol)
  }
})


# test_that("area of a square is calculated correctly", {
#   ball1 <- get_ball_mesh(m1)
#   m <- floor(m1 / 4)
#   square_l <- matrix(c(rep(-1, m), seq(-1, 1, length.out = m)),
#                      ncol = 2, byrow = FALSE)
#   square_r <- matrix(c(rep(1, m), seq(-1, 1, length.out = m)),
#                      ncol = 2, byrow = FALSE)
#   square_d <- matrix(c(seq(-1, 1, length.out = m), rep(-1, m)),
#                      ncol = 2, byrow = FALSE)
#   square_u <- matrix(c(seq(-1, 1, length.out = m), rep(1, m)),
#                      ncol = 2, byrow = FALSE)
#   square <- rbind(square_l, square_r, square_d, square_u)
#   err <- 4 - pi
#   err_est <- compute_error(square, ball1, m1, m2, function(x) 1)
#   print(paste0("Real: ", err))
#   print(paste0("Estimate: ", err_est))
#   expect_equal(err, err_est, tolerance = tol)
# })
# 
# 
# test_that("area of an ellipse is calculated correctly", {
#   ball1 <- get_ball_mesh(m1)
#   ball2 <- 10 * ball1
#   sigma1 <- 2 * matrix(c(11, 10.5, 10.5, 11.25), byrow = TRUE, ncol = 2)
#   lambda1 <- sqrtmat(sigma1)
#   ellipse1 <- ball1 %*% t(lambda1)
# 
#   sigma2 <- matrix(c(11, 10.5, 10.5, 11.25), byrow = TRUE, ncol = 2)
#   lambda2 <- sqrtmat(sigma2)
#   ellipse2 <- ball1 %*% t(lambda2)
# 
#   err <- pi * (prod(eigen(lambda1)$values) - prod(eigen(lambda2)$values))
#   err_est <- compute_error(ellipse1, ellipse2, m1, m2, function(x) 1)
#   print(paste0("Real: ", err))
#   print(paste0("Estimate: ", err_est))
#   expect_equal(err, err_est, tolerance = tol)
# 
#   err <- pi * (prod(eigen(lambda1)$values) - 1)
#   err_est <- compute_error(ellipse1, ball1, m1, m2, function(x) 1)
#   print(paste0("Real: ", err))
#   print(paste0("Estimate: ", err_est))
#   expect_equal(err, err_est, tolerance = tol)
# 
#   err <- pi * (10^2 - prod(eigen(lambda1)$values))
#   err_est <- compute_error(ellipse1, ball2, m1, m2, function(x) 1)
#   print(paste0("Real: ", err))
#   print(paste0("Estimate: ", err_est))
#   expect_equal(err, err_est, tolerance = tol)
# })
# 
# test_that("Error computation works for spherical t-distributions, when
#           estimate is far away from real", {
#   sigma <- diag(2)
#   p <- c(1 / 1000, 1 / 100)
#   err <- abs(p[1] - p[2])
# 
#   for (gamma_i in gamma) {
#     f <- function(x) dmvt(x, mu, sigma, 1 / gamma_i, log = FALSE)
#     real <- tdist_extreme_region(sigma, gamma_i, p[1], m1)
#     estimate <- tdist_extreme_region(sigma, gamma_i, p[2], m1)
#     err_est <- compute_error(real, estimate, m1, m2, f)
#     print(paste0("Real: ", err / p[1]))
#     print(paste0("Estimate: ", err_est / p[1]))
#     expect_equal(err / p[1], err_est / p[1], tolerance = tol)
#   }
# })
# 
# 
# test_that("Error computation works for elliptical t-distributions, when
#           estimate is far away from real", {
#   sigma <- matrix(c(11, 10.5, 10.5, 11.25), byrow = TRUE, ncol = 2)
#   p <- c(1 / 1000, 1 / 100)
#   err <- abs(p[1] - p[2])
# 
#   for (gamma_i in gamma) {
#     f <- function(x) dmvt(x, mu, sigma, 1 / gamma_i, log = FALSE)
#     real <- tdist_extreme_region(sigma, gamma_i, p[1], m1)
#     estimate <- tdist_extreme_region(sigma, gamma_i, p[2], m1)
#     err_est <- compute_error(real, estimate, m1, m2, f)
#     print(paste0("Real: ", err / p[1]))
#     print(paste0("Estimate: ", err_est / p[1]))
#     expect_equal(err / p[1], err_est / p[1], tolerance = tol)
#   }
# })
# 
# 
# test_that("Error computation works for spherical t-distributions when
#           estimate is close to real", {
#   sigma <- diag(2)
#   p <- c(1 / 1000, 1 / 1500)
#   err <- abs(p[1] - p[2])
# 
#   for (gamma_i in gamma) {
#     f <- function(x) dmvt(x, mu, sigma, 1 / gamma_i, log = FALSE)
#     real <- tdist_extreme_region(sigma, gamma_i, p[1], m1)
#     estimate <- tdist_extreme_region(sigma, gamma_i, p[2], m1)
#     err_est <- compute_error(real, estimate, m1, m2, f)
#     print(paste0("Real: ", err / p[1]))
#     print(paste0("Estimate: ", err_est / p[1]))
#     expect_equal(err / p[1], err_est / p[1], tolerance = tol)
#   }
# })
# 
# 
# test_that("Error computation works for elliptical t-distributions, when
#           estimate is close to real", {
#   sigma <- matrix(c(11, 10.5, 10.5, 11.25), byrow = TRUE, ncol = 2)
#   p <- c(1 / 1000, 1 / 1500)
#   err <- abs(p[1] - p[2])
# 
#   for (gamma_i in gamma) {
#     f <- function(x) dmvt(x, mu, sigma, 1 / gamma_i, log = FALSE)
#     real <- tdist_extreme_region(sigma, gamma_i, p[1], m1)
#     estimate <- tdist_extreme_region(sigma, gamma_i, p[2], m1)
#     err_est <- compute_error(real, estimate, m1, m2, f)
#     print(paste0("Real: ", err / p[1]))
#     print(paste0("Estimate: ", err_est / p[1]))
#     expect_equal(err / p[1], err_est / p[1], tolerance = tol)
#   }
# })
