library(extreme)
library(ggplot2)
library(ggforce)
library(latex2exp)
library(optparse)

#' Calculate angle between two vectors
#' @param a A vector
#' @param b A vector
#' @return Angle as radians
angle <- function(a, b) {
  acos(sum(a * b) / (norm(a, type = "2") * norm(b, type = "2")))
}

#' Calculate semiaxes lengths and rotation angle of ellipse of
#' form x^T \Sigma x = 1
#' @param mu a d-dimensional location of the ellipse
#' @param sigma inverted matrix corresponding to ellipse
#' @return vector of axis lengths, location and angle with respect to
#' positive x-axis
ellipse_to_vector <- function(mu, sigma) {
  # eigenvalues and eigenvectors
  decomp <- eigen(sigma)
  lambda <- decomp$values
  v <- decomp$vectors
  
  # lengths of the semiaxes
  semiaxes <- 1 / sqrt(lambda)
  # direction of major axis
  v_major  <- as.vector(v[, 2])
  
  # angle between major axis and x-axis
  theta <- angle(v_major, c(1, 0))
  if (v_major[2] < 0) {
    theta <- -theta
  }
  c(mu[1], mu[2], semiaxes[2], semiaxes[1], theta)
}


plot_ellipse <- function(tsample, q, qbar, qhat) {
  # Reshape data for ggplot
  q_v <- ellipse_to_vector(q$location, q$scale^(-2) * solve(q$scatter))
  qbar_v <- ellipse_to_vector(qbar$location,
                              qbar$scale^(-2) * solve(qbar$scatter))
  qhat_v <- ellipse_to_vector(qhat$location,
                              qhat$scale^(-2) * solve(qhat$scatter))
  data <- rbind(q_v, qbar_v, qhat_v, deparse.level = 0)
  data <- data.frame(x = data[, 1], y = data[, 2], a = data[, 3],
                     b = data[, 4], theta = data[, 5],
                     type = c("real", "bar", "hat"))
  tdf <- data.frame(x = tsample[, 1], y = tsample[, 2])
  
  # Plotting
  dir.create("figures", showWarnings = FALSE)
  f <- sprintf("%s/fig-n_%d-k_%d-p_%s-gamma_%.2f.eps", "figures", opt$n,
               opt$k, opt$p, opt$gamma)
  ggplot(data) +
    geom_ellipse(aes(x0 = x, y0 = y, a = a, b = b, angle = theta,
                     linetype = type), key_glyph = draw_key_path) +
    geom_point(data = tdf, aes(x = x, y = y)) +
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
                          values = c("real" = "solid", "hat" = "dashed",
                                     "bar" = "dotted"),
                          labels = c(TeX("$Q_p$"),
                                     TeX("$\\hat{Q}_p$"),
                                     TeX("$\\bar{Q}_p$")))
  ggsave(f, width = 7, height = 7)
} 


###############################################################################
#                                                                             #
#                           ARGUMENT PROCESSING                               #
#                                                                             #
###############################################################################

# Global constants
mu <- c(-25, 10)
sigma <- matrix(c(11, 10.5, 10.5, 11.25), byrow = TRUE, ncol = 2)
d <- 2

# Argument list
option_list <- list(
  make_option("--n", type = "integer", default = 5000,
              help = "Sample size"),
  make_option("--k", type = "numeric", default = 400,
              help = "Portion related to k"),
  make_option("--p", type = "character", default = "low",
              help = "Probability corresponding to quantile,
              choices are 'low', 'medium' and 'high'"),
  make_option("--gamma", type = "numeric", default = 1,
              help = "Choose extreme value index"),
  make_option("--seed", type = "integer", default = 204,
              help = "Set seed for sampling")
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Check arguments are valid
if (NA %in% opt) {
  stop("Argument is missing")
}

# Set parameter values
p <- switch(opt$p,
            low = 1 / 2000,
            medium = 1 / 5000,
            high = 1 / 10000,
            stop("Invalid value of 'p'."))
df <- 1 / opt$gamma

###############################################################################
#                                                                             #
#                           SIMULATION ROUTINE                                #
#                                                                             #
###############################################################################

# Generate sample
set.seed(opt$seed)
tsample <- relliptical(sqrt(d * rf(opt$n, d, df)), mu, sigma)

# Calculate real quantile region and estimates
q <- ellipsoidq(mu, sigma, sqrt(d * qf(1 - p, d, df)))
qbar <- qreg(tsample, p, method = "mcd", qmethod = "sample", alpha = 0.5)
qhat <- qreg(tsample, p, method = "mcd", qmethod = "extreme", k = opt$k,
             alpha = 0.5)

# Plotting
plot_ellipse(tsample, q, qbar, qhat)
