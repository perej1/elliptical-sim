# Estimate elliptical extreme quantile region for clover distribution with
# different parameters.
cli::cli_h1("Simulation starting")

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
  cli::cli_alert_success("Done")
}

cli::cli_h1("Simulation done")
