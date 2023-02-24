library(optparse)
library(rlang)
library(mvtnorm)
library(tibble)
library(readr)
library(pracma)

#' Generate Sample from Clover Distribution
#'
#' @param n An Integer. Sample size.
#'
#' @return A n by 2 matrix. Each row represents one observation.
gen_clover <- function(n) {
  u <- runif(n)
  v <- runif(n)
  r <- ((1 - u)^(-2) - 1)^(1 / 6)
  theta <- rep(0, n)
  for (i in 1:n) {
    f <- function(x) {
      (3 / (10 * pi) * r[i]^5 / (1 + r[i]^6)^(3 / 2) * (5 * x + sin(4 * x))) /
      (3 * r[i]^5 / (1 + r[i]^6)^(3 / 2)) - v[i]
    }
    theta[i] <- fzero(f, 0)$x
  }
  
  r * cbind(cos(theta), sin(theta))
}

# Argument list
option_list <- list(
  make_option("--type", type = "character", default = "cauchy",
              help = "Distribution type"),
  make_option("--n", type = "integer", default = 4,
              help = "Sample size"),
  make_option("--d", type = "integer", default = 2,
              help = "Dimension of distribution"),
  make_option("--s", type = "integer", default = 5,
              help = "Number of samples"),
  make_option("--seed", type = "integer", default = 204,
              help = "Set seed for sampling")
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Nontrivial location-scatter pair
mu <- c(100, -250)
sigma <- matrix(c(11, 10.5, 10.5, 11.25), byrow = TRUE, ncol = 2)

# Simulate s samples from distribution specified by type
gen_samp <- function(n) {
  switch (opt$type,
    cauchy = rmvt(n, diag(2), 1, rep(0, 2)),
    cauchyAff = rmvt(n, sigma, 1, mu),
    cauchy3d = rmvt(n, diag(3), 1, rep(0, 3)),
    tdistDeg2 = rmvt(n, diag(2), 2, rep(0, 2)),
    tdistDeg4 = rmvt(n, diag(2), 4, rep(0, 2)),
    clover = gen_clover(n),
    abort("Invalid distribution type.")
  )
}

set.seed(opt$seed)
samples <- replicate(opt$s, gen_samp(opt$n), simplify = FALSE)
samples <- do.call(cbind, samples)

if (opt$d == 2) {
  colnames(samples) <- paste0(c("x", "y", "z"), rep(1:opt$s, each = 3))
} else if (opt$d == 3) {
  colnames(samples) <- paste0(c("x", "y"), rep(1:opt$s, each = 2))
} else {
  abort("Dimension must be two or three.")
}

# Write data
samples <- as_tibble(samples)
filename <- paste0(opt$type, "_", opt$n, ".csv")
write_csv(samples, paste0("data/samples/", filename))
