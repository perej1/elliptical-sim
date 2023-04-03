# Summarise simulation results for a particular scenario
source("functions.R")

option_list <- list(
  make_option("--type", type = "character", default = "cauchyAff",
              help = "Distribution type"),
  make_option("--n", type = "integer", default = 500,
              help = "Sample size"),
  make_option("--p", type = "character", default = "low",
              help = "Probability mass outside quantile region"),
  make_option("--k", type = "character", default = "small",
              help = "Sample size of the tail"),
  make_option("--seed", type = "integer", default = 278,
              help = "Set seed for sampling")
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Read data
filename <- paste0("type_", opt$type,
                   "_n_", opt$n,
                   "_seed_", opt$seed,
                   "_p_", opt$p,
                   "_k_", opt$k, ".csv"
                   )
elliptical_estimates <- readr::read_csv(paste0("sim-data/elliptical-estimates/",
                                               filename),
                                        show_col_types = FALSE)
depth_estimates <- readr::read_csv(paste0("sim-data/depth-estimates/",
                                          filename),
                                   show_col_types = FALSE)
errors <- readr::read_csv(paste0("sim-data/errors/",
                                 filename),
                          show_col_types = FALSE)

filename <- paste0("type_", opt$type, "_p_", opt$p, ".csv")
real <- readr::read_csv(paste0("sim-data/real-regions/", filename),
                        show_col_types = FALSE)

# Construct table of median, max and min errors
mins <- apply(errors, 2, min, na.rm = TRUE)
maxs <- apply(errors, 2, max, na.rm = TRUE)
medians <- apply(errors, 2, median, na.rm = TRUE)

stats_df <- tibble::as_tibble(rbind(mins, maxs, medians)) %>%
  mutate(type = c("min", "max", "median"))

# Extract estimates corresponding to min, max and median errors.
ellipse_min <- elliptical_estimates %>%
  select(num_range(c("x", "y"), which.min(errors$elliptical))) %>%
  rename(x = 1, y = 2)
ellipse_max <- elliptical_estimates %>%
  select(num_range(c("x", "y"), which.max(errors$elliptical))) %>%
  rename(x = 1, y = 2)
ellipse_median <- elliptical_estimates %>%
  select(num_range(c("x", "y"), 
                   which.min(abs(errors$elliptical - medians["elliptical"])))
         ) %>%
  rename(x = 1, y = 2)

depth_min <- depth_estimates %>%
  select(num_range(c("x", "y"), which.min(errors$depth))) %>%
  rename(x = 1, y = 2)
depth_max <- depth_estimates %>%
  select(num_range(c("x", "y"), which.max(errors$depth))) %>%
  rename(x = 1, y = 2)
depth_median <- depth_estimates %>%
  select(num_range(c("x", "y"), 
                   which.min(abs(errors$depth - medians["depth"])))
  ) %>%
  rename(x = 1, y = 2)

# Combine data for plotting
data_min <- bind_rows(real, ellipse_min, depth_min) %>%
  mutate(group = rep(c("real", "ellipse", "depth"), each = 1000))

data_max <- bind_rows(real, ellipse_max, depth_max) %>%
  mutate(group = rep(c("real", "ellipse", "depth"), each = 1000))

data_median <- bind_rows(real, ellipse_median, depth_median) %>%
  mutate(group = rep(c("real", "ellipse", "depth"), each = 1000))

# Plotting
g_min <- plot_data(data_min)
g_max <- plot_data(data_max)
g_median <- plot_data(data_median)

# Save results
filename <- paste0("type_", opt$type,
                   "_n_", opt$n,
                   "_seed_", opt$seed,
                   "_p_", opt$p,
                   "_k_", opt$k)
readr::write_csv(stats_df,
                 paste0("summary-data/stats/", filename, ".csv")
)
ggsave(paste0("summary-data/figures-min/", filename, ".jpg"),
       plot = g_min, width = 7, height = 7, dpi = 1000)
ggsave(paste0("summary-data/figures-max/", filename, ".jpg"),
       plot = g_max, width = 7, height = 7, dpi = 1000)
ggsave(paste0("summary-data/figures-median/", filename, ".jpg"),
       plot = g_median, width = 7, height = 7, dpi = 1000)
