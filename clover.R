m <- 1000
beta <- 0.01

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

clover_contour <- function(beta, m) {
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
    r[i] <- pracma::fzero(g, c(0, 500))$x
  }
  
  x <- r * cos(theta)
  y <- r * sin(theta)
  cbind(x, y)
}

plot(clover_contour(beta, m), type = "l")
