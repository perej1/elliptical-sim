library(optparse)
library(ggplot2)
suppressPackageStartupMessages(library(dplyr))


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


elliptical_extreme_qregion <- function(data, sigma_est, p, k, m) {
  n <- nrow(data)
  w <- get_ball_mesh(m)
  
  # Approximate generating variate
  radius <- sqrt(stats::mahalanobis(data, FALSE, sigma_est, inverted = FALSE))
  radius_sort <- sort(radius, decreasing = FALSE)
  
  # Estimate extreme value index
  gamma_est <- mean((log(radius_sort[(n - k):n]) - log(radius_sort[n - k]))[-1])
  
  # Estimate extreme quantile of generating variate
  r_hat <- radius_sort[n - k] * (k / (n * p))^gamma_est
  
  # Estimate extreme quantile region
  lambda <- sqrtmat(r_hat^2 * sigma_est)
  w %*% t(lambda)
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


option_list <- list(
  make_option("--n", type = "integer", default = 5000,
              help = "Sample size"),
  make_option("--k", type = "character", default = "large",
              help = "Sample size of the tail")
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Global parameters
m <- 1000

# Set value for k and p
k <- switch(opt$k,
            large = 0.2 * opt$n,
            medium = 0.1 * opt$n,
            small = 0.05 * opt$n,
            rlang::abort("Invalid value of k")
)

p <- c(2 / opt$n, 1 / (2 * opt$n))

# Generate sample
set.seed(204)
sample <- gen_clover(opt$n)
colnames(sample) <- c("x", "y")

# Estimate scatter
est <- robustbase::covMcd(sample, alpha = 0.5)

data_list <- as.list(rep(NA, length(p)))
for (i in 1:length(p)) {
  estimate <- elliptical_extreme_qregion(sample, est$cov, p[i], k, m)
  real <- clover_contour_p(p[i], m)
  
  colnames(estimate) <- c("x", "y")
  colnames(real) <- c("x", "y")
  
  data_list[[i]] <- bind_rows(as_tibble(real), as_tibble(estimate)) %>%
    mutate(group = rep(paste0(c("real", "estimate"), i), each = 1000))
}


g<-ggplot(as_tibble(sample), aes(x=x, y=y)) +
  geom_point() +
  geom_path(data = bind_rows(data_list),
            aes(x=x, y=y, group = group, linetype = group),
            show.legend = FALSE) +
  coord_fixed() +
  theme(axis.title = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.key.size = unit(1, "cm"),
        legend.text = element_text(size = 15),
        axis.text = element_text(size = 15)) +
  scale_linetype_manual(name = NULL,
                        values = c("real1" = "solid", "estimate1" = "dashed",
                                   "real2" = "solid", "estimate2" = "dashed"))

# Save figure
filename <- paste0("n_", opt$n, "_k_", opt$k, ".jpg")
ggsave(paste0("clover-figures/", filename), plot = g, width = 7, height = 7,
       dpi = 1000)
