library(ggplot2)
library(optparse)
suppressPackageStartupMessages(library(dplyr))

# Argument list
option_list <- list(
  make_option("--type", type = "character", default = "cauchy",
              help = "Distribution type"),
  make_option("--n", type = "integer", default = 1000,
              help = "Sample size"),
  make_option("--p", type = "character", default = "low",
              help = "Probability mass outside quantile region"),
  make_option("--k", type = "character", default = "large",
              help = "Sample size of the tail"),
  make_option("--i", type = "integer", default = 2,
              help = "Scenario number")
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
m <- 1000

# Read data about estimates and theoretical quantile region
filename <- paste0(opt$type, "_", opt$n, "_", opt$p, "_", opt$k, ".csv")

ellipse_i <- readr::read_csv(paste0("data/elliptical-estimates/", filename),
                             show_col_types = FALSE) %>%
  select(num_range(c("x", "y"), opt$i)) %>%
  rename(x = 1, y = 2)

depth_i <- readr::read_csv(paste0("data/depth-estimates/", filename),
                           show_col_types = FALSE) %>%
  select(num_range(c("x", "y"), opt$i)) %>%
  rename(x = 1, y = 2)

filename <- paste0(opt$type, "_", opt$n, "_", opt$p, ".csv")
real <- readr::read_csv(paste0("data/real-regions/", filename),
                        show_col_types = FALSE)

data_est <- bind_rows(real, ellipse_i, depth_i) %>%
  mutate(group = rep(c("real", "ellipse", "depth"), each = 1000))

# Read data about sample
filename <- paste0(opt$type, "_", opt$n, ".csv")
sample_i <- readr::read_csv(paste0("data/samples/", filename),
                            show_col_types = FALSE) %>%
  select(num_range(c("x", "y"), opt$i)) %>%
  rename(x = 1, y = 2)


g <- ggplot(data_est, aes(x = x, y = y)) +
  geom_path(aes(group = group, linetype = group), show.legend = FALSE) +
  geom_point(data = sample_i, mapping = aes(x = x, y = y)) +
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
                        values = c("real" = "solid", "ellipse" = "dashed",
                                   "depth" = "dotted"))

# Save figure
filename <- paste0(opt$type, "_", opt$n, "_", opt$p, "_", opt$k, "_", opt$i, ".jpg")
ggsave(paste0("figures/", filename), plot = g, width = 7, height = 7, dpi = 1000)
