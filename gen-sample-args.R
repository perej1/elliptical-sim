# Make a txt file of simulation arguments, one argument combination per row
type <- c("cauchy", "cauchyAff", "cauchy3d", "tdistDeg4", "tdistDeg2")
n <- c(100, 500, 1000, 5000)
s <- 100
seed <- 204

arg <- expand.grid(type = type, n = n, s = s, seed = seed)
arg <- sprintf("gen-sample.R --type %s --n %d --s %d --seed %d",
               arg$type,
               arg$n,
               arg$s,
               arg$seed)

connection <- file("args/arg-gen-sample.txt")
writeLines(arg, connection)
close(connection)
