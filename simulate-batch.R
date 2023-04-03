# Perform simulations for different scenarios.
library(optparse)

option_list <- list(
  make_option("--simulate", type = "logical", default = FALSE,
              help = "Simulate or not"),
  make_option("--summarise", type = "logical", default = FALSE,
              help = "Summarise or not"),
  make_option("--clover", type = "logical", default = TRUE,
              help = "Generate clover example or not")
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

type <- c("cauchy", "cauchyAff", "tdistDeg4")
n <- c(1000, 5000)
p <- c("low", "medium", "high")
k <- c("small", "medium", "large")
seed <- 278

arg <- expand.grid(type = type, n = n, p = p, k = k, seed = seed)

if (opt$simulate) {
  cli::cli_h1("Simulation starting")
  arg_vector <- sprintf(paste0("Rscript simulate.R --type %s --n %d --p %s ",
                               "--k %s --seed %d"),
                        arg$type,
                        arg$n,
                        arg$p,
                        arg$k,
                        arg$seed)

  n_scenario <- length(arg_vector)
  for (i in 1:n_scenario) {
    cli::cli_h2("Scenario {i}/{n_scenario} running")
    cli::cli_alert_info("Parameters:
                      . \t type={arg$type[i]}
                      . \t n={arg$n[i]}
                      . \t p={arg$p[i]}
                      . \t k={arg$k[i]}
                      . \t seed={arg$seed[i]}")
    system(arg_vector[i])
  }

  cli::cli_h1("Simulation done")
}

if (opt$summarise) {
  cli::cli_h1("Summarising starting")
  arg_vector <- sprintf(paste0("Rscript summarise.R --type %s --n %d --p %s ",
                               "--k %s --seed %d"),
                        arg$type,
                        arg$n,
                        arg$p,
                        arg$k,
                        arg$seed)
  n_scenario <- length(arg_vector)
  cli::cli_progress_bar("Summarizing:", total = n_scenario)
  for (i in 1:n_scenario) {
    system(arg_vector[i])
    cli::cli_progress_update()
  }

  cli::cli_h1("Summary done")
}

if (opt$clover) {
  cli::cli_h1("Clover example starting")

  n <- 5000
  k <- c("small", "medium", "large")
  seed <- 278

  arg <- expand.grid(n = n, k = k, seed = seed)
  arg_vector <- sprintf("Rscript estimate-clover.R --n %d --k %s --seed %d",
                        arg$n,
                        arg$k,
                        arg$seed)

  n_scenario <- length(arg_vector)
  for (i in 1:n_scenario) {
    cli::cli_h2("Scenario {i}/{n_scenario} running")
    cli::cli_alert_info("Parameters:
                      . \t n={arg$n[i]}
                      . \t k={arg$k[i]}
                      . \t seed={arg$seed[i]}")
    system(arg_vector[i])
  }

  cli::cli_h1("Clover example done")
}
