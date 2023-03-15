# Perform simulations for different scenarios.
cli::cli_h1("Simulation starting")

type <- c("cauchy", "cauchyAff", "tdistDeg4", "tdistDeg2")
n <- c(500, 1000, 5000)
p <- c("low", "medium", "high")
k <- c("small", "medium", "large")
seed <- 278

arg <- expand.grid(type = type, n = n, p = p, k = k, seed = seed)
arg_vector <- sprintf("Rscript simulate.R --type %s --n %d --p %s --k %s --seed %d",
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
  cli::cli_alert_success("Done")
}

cli::cli_h1("Simulation done")
