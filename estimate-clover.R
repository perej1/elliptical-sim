# Estimate elliptical extreme quantile region for clover distribution with
# certain parameters. Save the resulting figure.
source("functions.R")

option_list <- list(
  make_option("--n", type = "integer", default = 5000,
              help = "Sample size"),
  make_option("--k", type = "character", default = "large",
              help = "Sample size of the tail"),
  make_option("--seed", type = "integer", default = 204,
              help = "Set seed for sampling")
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Global parameters
m <- 1000

# Set values for k and p
k <- switch(opt$k,
            large = 0.2 * opt$n,
            medium = 0.1 * opt$n,
            small = 0.05 * opt$n,
            rlang::abort("Invalid value of k")
)

p <- c(2 / opt$n, 1 / (2 * opt$n))

# Generate sample
set.seed(opt$seed)
sample <- gen_clover(opt$n)
colnames(sample) <- c("x", "y")

# Estimate scatter and location
est <- robustbase::covMcd(sample, alpha = 0.5)

# Compute real and estimated quantile regions
data_list <- as.list(rep(NA, length(p)))
for (i in 1:length(p)) {
  estimate <- elliptical_extreme_qregion(sample, est$center, est$cov, p[i], k, m)
  real <- clover_contour_p(p[i], m)
  
  colnames(estimate) <- c("x", "y")
  colnames(real) <- c("x", "y")
  
  data_list[[i]] <- bind_rows(as_tibble(real), as_tibble(estimate)) %>%
    mutate(group = rep(paste0(c("real", "estimate"), i), each = 1000))
}

# Plotting
g <- ggplot(as_tibble(sample), aes(x = x, y = y)) +
  geom_point() +
  geom_path(data = bind_rows(data_list),
            aes(x = x, y = y, group = group, linetype = group),
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
