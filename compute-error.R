library(optparse)
library(mvtnorm)

get_ball_mesh <- function(m) {
  w <- seq(0, 2 * pi, length.out = m)
  cbind(cos(w), sin(w))
}

sqrtmat <- function(sigma) {
  eigenval <- eigen(sigma)$values
  if (any(eigenval <= 0) || any(sigma != t(sigma))) {
    rlang::abort("`sigma` must be a symmetric positive definite matrix.")
  }
  eigenvec <- eigen(sigma)$vectors
  eigenvec %*% diag(eigenval^0.5) %*% t(eigenvec)
}

tdist_extreme_region <- function(mu, sigma, gamma, p, m) {
  w <- get_ball_mesh(m)
  lambda <- sqrtmat(2 * stats::qf(1 - p, 2, 1 / gamma) * sigma)
  sweep(w %*% t(lambda), 2, mu, "+")
}

ellipse_region <- function(mu, sigma, m) {
  w <- get_ball_mesh(m)
  lambda <- sqrtmat(sigma)
  sweep(w %*% t(lambda), 2, mu, "+")
}

compute_error <- function(real, estimate, m1, m2, f) {
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

###################### TEST ESTIMATION OF AREA ####################

m <- 10000
mu <- rep(0, 2)
sigma1 <- matrix(c(9, 1.5, 1.5, 5), byrow = TRUE, ncol = 2)
#sigma2 <- matrix(c(14, 0.5, 0.5, 2), byrow = TRUE, ncol = 2)
sigma2 <- matrix(c(11, 2, 2, 11.25), byrow = TRUE, ncol = 2)

e1 <- ellipse_region(mu, sigma1, m)
e2 <- ellipse_region(mu, sigma2, m)

comp <- compute_error(e1, e2, m, 1000, function (x) 1)
#xy <- comp$coord

plot(e2[, 1], e2[, 2], type = "l", asp = 1)
points(e1[, 1], e1[, 2], type = "l")
#points(xy[, 1], xy[, 2], pch = 21)

pi * abs(prod(eigen(sigma1)$values^(0.5)) - prod(eigen(sigma2)$values^(0.5)))
comp$res # Estimate

#################### TEST ESTIMATION OF DENSITY ######################

m <- 1000
#sigma <- diag(2)
#sigma <- matrix(c(11, 10.5, 10.5, 11.25), byrow = TRUE, ncol = 2)
sigma <- matrix(c(2, 0, 0, 1), byrow = TRUE, ncol = 2)
sigma <- matrix(c(11, 2, 2, 11.25), byrow = TRUE, ncol = 2)
#sigma <- matrix(c(4, 0, 0, 2), byrow = TRUE, ncol = 2)
gamma <- 1
p <- c(1 / 1000, 1 / 1500)

e1 <- tdist_extreme_region(mu, sigma, gamma, p[1], m)
e2 <- tdist_extreme_region(mu, sigma, gamma, p[2], m)

# Real prob mass
p_diff_real <- (1 - p[2]) - (1 - p[1])

f <- function(x) {
  dmvt(x, mu, sigma, df = 1 / gamma, log = FALSE)
}

comp <- compute_error(e1, e2, m, 100, f)
xy <- comp$coord

plot(e2[, 1], e2[, 2], type = "l", asp = 1)
points(e1[, 1], e1[, 2], type = "l")
points(xy[, 1], xy[, 2], pch = 21)

p_diff_real
comp$res
