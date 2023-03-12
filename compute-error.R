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

compute_error <- function(r, r_hat, theta, m2, f) {
  m1 <- length(theta)
  res <- 0
  coord <- matrix(NA, ncol = 2, nrow = m1 * m2)
  
  for (i in 1:m1) {
    r_seq <- seq(min(r_hat[i], r[i]), max(r_hat[i], r[i]), length.out = m2)
    x <- r_seq * cos(theta[i])
    y <- r_seq * sin(theta[i])
    coord[((i-1)*m2+1):(i*m2),] <- cbind(x,y)
    resi <- 0
    for (j in 1:m2) {
      resi <- resi + f(c(x[j], y[j])) * r_seq[j]
    }
    res <- res + abs(r[i] - r_hat[i]) * resi
  }
  list(res = 2 * pi / (m1 * m2) * res, coord = coord)
}

###################### TEST AREA ####################

m <- 1000
mu <- rep(0, 2)
sigma1 <- matrix(c(1, 0, 0, 1), byrow = TRUE, ncol = 2)
sigma2 <- matrix(c(4, 0, 0, 1), byrow = TRUE, ncol = 2)

e1 <- ellipse_region(mu, sigma1, m)
e2 <- ellipse_region(mu, sigma2, m)

plot(e2[, 1], e2[, 2], type = "l", asp = 1)
points(e1[, 1], e1[, 2], type = "l")

# Real area
a_real <- pi * 2 - pi

# Radius
r_e1 <- apply(e1, 1, norm, type = "2")
r_e2 <- apply(e2, 1, norm, type = "2")

# Angle
theta_e1 <- rep(NA, m)
theta_e2 <- rep(NA, m)

e11 <- sweep(e1, 1, r_e1, "/")
e22 <- sweep(e2, 1, r_e2, "/")

for (i in 1:m) {
  #theta_e1[i] <- acos(e1[i, 1] / r_e1[i])
  #if (e1[i, 2] < 0) {
  #  theta_e1[i] <- 2 * pi - theta_e1[i]
  #}
  
  theta_e1[i] <- acos(e11[i, 1])
  if (e11[i, 2] < 0) {
    theta_e1[i] <- 2 * pi - theta_e1[i]
  }
  
  # theta_e2[i] <- acos(e2[i, 1] / r_e2[i])
  # if (e2[i, 2] < 0) {
  #   theta_e2[i] <- 2 * pi - theta_e2[i]
  # }
  
  theta_e2[i] <- acos(e22[i, 1])
  if (e22[i, 2] < 0) {
    theta_e2[i] <- 2 * pi - theta_e2[i]
  }
}

e1_df <- data.frame(x = e1[, 1], y = e1[, 2], r = r_e1, theta = theta_e1)
e2_df <- data.frame(x = e2[, 1], y = e2[, 2], r = r_e2, theta = theta_e2)

ind <- apply(e11 %*% t(e22), 1, which.max)
#ind2 <- apply(e1 %*% t(e2), 1, which.max)

e2_df <- e2_df[ind, ]

# Test that angle is calculated correctly
test1 <- e1_df$r * cbind(cos(e1_df$theta), sin(e1_df$theta))
test2 <- e2_df$r * cbind(cos(e1_df$theta), sin(e1_df$theta))

plot(test2[, 1], test2[, 2], type = "l", asp = 1)
points(test1[, 1], test1[, 2], type = "l")

plot(e1_df$x, e1_df$y, type = "l", asp = 1)
points(e2_df$x, e2_df$y, type = "l")

plot(e2_df$x, e2_df$y, type = "l", asp = 1)
points(test2[, 1], test2[, 2], type = "l")

f <- function(x) {
  1
}

comp <- compute_error(e1_df$r, e2_df$r, e1_df$theta, 40, f)
xy <- comp$coord

plot(test2[, 1], test2[, 2], type = "l", asp = 1)
points(test1[, 1], test1[, 2], type = "l")
points(xy[,1], xy[,2], pch = 21)

a_real
comp$res

#################### TEST DENS ######################
# Test the function
m <- 100
mu <- rep(0, 2)
sigma <- diag(2)
sigma <- matrix(c(11, 10.5, 10.5, 11.25), byrow = TRUE, ncol = 2)
gamma <- 1
p <- c(1 / 1000, 1 / 1500)

e1 <- tdist_extreme_region(mu, sigma, gamma, p[1], m)
e2 <- tdist_extreme_region(mu, sigma, gamma, p[2], m)

ind <- sample(1:m, 10, replace = FALSE)
plot(e2[, 1], e2[, 2], type = "l", asp = 1)
points(e1[, 1], e1[, 2], type = "l")

points(e2[ind, 1], e2[ind, 2], pch = 21)
points(e1[ind, 1], e1[ind, 2], pch = 21)

# Real prob mass
p_diff_real <- (1 - p[2]) - (1 - p[1])

# Radius
r_e1 <- apply(e1, 1, norm, type = "2")
r_e2 <- apply(e2, 1, norm, type = "2")

# Angle
theta_e1 <- rep(NA, m)
theta_e2 <- rep(NA, m)

e11 <- sweep(e1, 1, r_e1, "/")
e22 <- sweep(e2, 1, r_e2, "/")

for (i in 1:m) {
  #theta_e1[i] <- acos(e1[i, 1] / r_e1[i])
  #if (e1[i, 2] < 0) {
  #  theta_e1[i] <- 2 * pi - theta_e1[i]
  #}
  
  theta_e1[i] <- acos(e11[i, 1])
  if (e11[i, 2] < 0) {
   theta_e1[i] <- 2 * pi - theta_e1[i]
  }
  
  theta_e2[i] <- acos(e2[i, 1] / r_e2[i])
  if (e2[i, 2] < 0) {
    theta_e2[i] <- 2 * pi - theta_e2[i]
  }
}

test1 <- r_e1 * cbind(cos(theta_e2), sin(theta_e2))
test2 <- r_e2 * cbind(cos(theta_e2), sin(theta_e2))
plot(test2[, 1], test2[, 2], type = "l", asp = 1)
points(test1[, 1], test1[, 2], type = "l")

f <- function(x) {
  dmvt(x, mu, sigma, df = 1 / gamma, log = FALSE)
}

comp <- compute_error(r_e1, r_e2, theta_e1, 40, f)
xy <- comp$coord

plot(test2[, 1], test2[, 2], type = "l", asp = 1)
points(test1[, 1], test1[, 2], type = "l")
points(xy[,1], xy[,2], pch = 21)

p_diff_real
comp$res
