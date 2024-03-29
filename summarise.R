source("functions.R")

# Read argument combinations
args <- readr::read_csv("sim-args.csv", show_col_types = FALSE)

args_2d <- args %>%
  filter(d == 2 & type != "tdistSkew")

args_3d <- args %>%
  filter(d == 3)

medians_3d <- args_3d %>%
  mutate(elliptical = NA, depth = NA)

args_skew <- args %>%
  filter(type == "tdistSkew")

# Compute min, max and median errors for 2d and 3d cases.
# Produce corresponding figures for the 2d case.
# Save medians for the 3d cases.
for (i in 1:(nrow(args_2d) + nrow(args_3d))) {
  arg <- rbind(args_3d, args_2d)[i, ]
  filename <- paste0("type_", arg$type,
                     "_s_", arg$s,
                     "_d_", arg$d,
                     "_n_", arg$n,
                     "_p_", arg$p,
                     "_k_", arg$k,
                     "_seed_", arg$seed)
  errors <- readr::read_csv(paste0("sim-data/errors/", filename, ".csv"),
                            show_col_types = FALSE)
  
  mins <- apply(errors, 2, min, na.rm = TRUE)
  maxs <- apply(errors, 2, max, na.rm = TRUE)
  medians <- apply(errors, 2, median, na.rm = TRUE)
  tibble::as_tibble(rbind(mins, maxs, medians)) %>%
    mutate(type = c("min", "max", "median")) %>%
    readr::write_csv(paste0("summary-data/stats/", filename, ".csv"))
  
  if (arg$d == 3) {
    medians_3d[i, c("elliptical", "depth")] <- list(medians[1], medians[2])
  }
  
  if (arg$d == 2) {
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

# Plot scenarios for skewed t distribution.
# Each plot includes cases for all p.
args_skew_no_p <- args_skew %>%
  select(-p) %>%
  distinct()
p <- c("low", "high")

for (i in 1:nrow(args_skew_no_p)) {
  arg <- args_skew_no_p[i, ]
  estimate_list <- vector("list", length(p))
  real_list <- vector("list", length(p))
  
  filename_sample <- paste0("type_", arg$type,
                            "_s_", arg$s,
                            "_d_", arg$d,
                            "_n_", arg$n,
                            "_seed_", arg$seed, ".csv")
  sample <- readr::read_csv(paste0("sim-data/samples/", filename_sample),
                            show_col_types = FALSE) %>%
    rename(x = x1, y = y1)
  
  for (j in seq_along(p)) {
    filename_real <- paste0("type_", arg$type,
                            "_d_", arg$d,
                            "_n_", arg$n,
                            "_p_", p[j], ".csv")
    real <- readr::read_csv(paste0("sim-data/real-regions/", filename_real),
                            show_col_types = FALSE)
    real_list[[j]] <- real
    
    filename_estimate <- paste0("type_", arg$type,
                                "_s_", arg$s,
                                "_d_", arg$d,
                                "_n_", arg$n,
                                "_p_", p[j],
                                "_k_", arg$k,
                                "_seed_", arg$seed, ".csv")
    estimate <- readr::read_csv(paste0("sim-data/elliptical-estimates/",
                                       filename_estimate),
                                show_col_types = FALSE) %>%
      rename(x = x1, y = y1)
    estimate_list[[j]] <- estimate
  }
  
  g <- plot_sample_and_estimates(sample, real_list, estimate_list, p)
  filename <- paste0("type_", arg$type,
                     "_s_", arg$s,
                     "_d_", arg$d,
                     "_n_", arg$n,
                     "_k_", arg$k,
                     "_seed_", arg$seed)
  ggsave(paste0("summary-data/figures-skew/", filename, ".jpg"),
         plot = g, width = 7, height = 7, dpi = 1000)
}


# Plot scenarios for 3d cases.
# Each plot includes boxplots for errors of depth/elliptical estimators for
# optimal k with certain combination of n and p.
n <- c(1000, 5000)
p <- c("low", "medium", "high")
args_3d_no_pn <- args_3d %>%
  select(-p, -n) %>%
  distinct()

for (i in 1:nrow(args_3d_no_pn)) {
  arg <- args_3d_no_pn[i, ]
  for (j1 in seq_along(p)) {
    for (j2 in seq_along(n)) {
      medians <- filter(medians_3d, type == arg$type & p == p[j1] & n == n[j2])
      
      k_elliptical <- medians$k[which.min(medians$elliptical)]
      filename_elliptical <- paste0("type_", arg$type,
                                    "_s_", arg$s,
                                    "_d_", arg$d,
                                    "_n_", n[j2],
                                    "_p_", p[j1],
                                    "_k_", k_elliptical,
                                    "_seed_", arg$seed, ".csv")
      elliptical <- readr::read_csv(paste0("sim-data/errors/",
                                           filename_elliptical),
                                    show_col_types = FALSE) %>%
        select(elliptical)
      
      k_depth <- medians$k[which.min(medians$depth)]
      filename_depth <- paste0("type_", arg$type,
                                    "_s_", arg$s,
                                    "_d_", arg$d,
                                    "_n_", n[j2],
                                    "_p_", p[j1],
                                    "_k_", k_depth,
                                    "_seed_", arg$seed, ".csv")
      depth <- readr::read_csv(paste0("sim-data/errors/", filename_elliptical),
                               show_col_types = FALSE) %>%
        select(depth)
      
      g <- boxplot_errors(elliptical, depth)
      
      filename <- paste0("type_", arg$type,
                         "_s_", arg$s,
                         "_d_", arg$d,
                         "_n_", n[j2],
                         "_p_", p[j1],
                         "_seed_", arg$seed)
      ggsave(paste0("summary-data/figures-boxplot/", filename, ".jpg"),
             plot = g, width = 7, height = 7, dpi = 1000)
    }
  }
}
