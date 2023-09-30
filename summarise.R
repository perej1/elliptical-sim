source("functions.R")

# Read argument combinations
args <- readr::read_csv("sim-args.csv", show_col_types = FALSE)

args_2d <- args %>%
  filter(d == 2 & type != "tdistSkew")

args_3d <- args %>%
  filter(d == 3)

args_skew <- args %>%
  filter(type == "tdistSkew")

# Compute min, max and median errors for 2d and 3d cases.
# Produce corresponding figures for the 2d case.
for (i in 1:(nrow(args_2d) + nrow(args_3d))) {
  # Read error data
  arg <- rbind(args_2d, args_3d)[i, ]
  filename <- paste0("type_", arg$type,
                     "_s_", arg$s,
                     "_d_", arg$d,
                     "_n_", arg$n,
                     "_p_", arg$p,
                     "_k_", arg$k,
                     "_seed_", arg$seed)
  errors <- readr::read_csv(paste0("sim-data/errors/", filename, ".csv"),
                            show_col_types = FALSE)
  
  # Construct table of median, max and min errors and write data
  mins <- apply(errors, 2, min, na.rm = TRUE)
  maxs <- apply(errors, 2, max, na.rm = TRUE)
  medians <- apply(errors, 2, median, na.rm = TRUE)
  tibble::as_tibble(rbind(mins, maxs, medians)) %>%
    mutate(type = c("min", "max", "median")) %>%
    readr::write_csv(paste0("summary-data/stats/", filename, ".csv"))
  
  if (arg$d == 2) {
    # Read data corresponding estimates and theoretical quantile regions
    elliptical_estimates <- readr::read_csv(paste0("sim-data/elliptical-estimates/",
                                                   filename, ".csv"),
                                            show_col_types = FALSE)
    depth_estimates <- readr::read_csv(paste0("sim-data/depth-estimates/",
                                              filename, ".csv"),
                                       show_col_types = FALSE)
    
    filename_real <- paste0("type_", arg$type,
                            "_d_", arg$d,
                            "_n_", arg$n,
                            "_p_", arg$p, ".csv")
    real <- readr::read_csv(paste0("sim-data/real-regions/", filename_real),
                            show_col_types = FALSE)
    
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
    # data_min <- bind_rows(real, ellipse_min, depth_min) %>%
    #   mutate(group = rep(c("real", "ellipse", "depth"), each = 100))
    # 
    # data_max <- bind_rows(real, ellipse_max, depth_max) %>%
    #   mutate(group = rep(c("real", "ellipse", "depth"), each = 100))
    # 
    # data_median <- bind_rows(real, ellipse_median, depth_median) %>%
    #   mutate(group = rep(c("real", "ellipse", "depth"), each = 100))
    
    # Plotting
    g_min <- plot_real_estimate(real, ellipse_min, depth_min)
    g_max <- plot_real_estimate(real, ellipse_max, depth_max)
    g_median <- plot_real_estimate(real, ellipse_median, depth_median)
    
    ggsave(paste0("summary-data/figures-min/", filename, ".jpg"),
           plot = g_min, width = 7, height = 7, dpi = 1000)
    ggsave(paste0("summary-data/figures-max/", filename, ".jpg"),
           plot = g_max, width = 7, height = 7, dpi = 1000)
    ggsave(paste0("summary-data/figures-median/", filename, ".jpg"),
           plot = g_median, width = 7, height = 7, dpi = 1000)
  }
}
