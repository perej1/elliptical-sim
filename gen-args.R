# Simulation arguments
type <- c("cauchy", "cauchyAff")
#type <- c("cauchy", "cauchyAff", "cauchy3d", "tdistDeg4", "tdistDeg2")
n <- c(100, 500, 1000, 5000)
seed <- 204
p <- c("low", "medium", "high")
k <-c("large", "medium", "small")

# Simulation arguments for creating samples
arg <- expand.grid(type = type, n = n, seed = seed)
arg <- sprintf("gen-sample.R --type %s --n %d --seed %d",
               arg$type,
               arg$n,
               arg$seed)

connection <- file("args/arg-gen-sample.txt")
writeLines(arg, connection)
close(connection)

# Simulation arguments for estimation
arg <- expand.grid(type = type, n = n, p = p, k = k)
arg <- sprintf("estimate.R --type %s --n %d --p %s --k %s",
               arg$type,
               arg$n,
               arg$p,
               arg$k)

connection <- file("args/arg-estimate.txt")
writeLines(arg, connection)
close(connection)
