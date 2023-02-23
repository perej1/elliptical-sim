# Make a txt file of simulation arguments, one argument combination per row
n <- c(100, 500, 1000, 5000)
k <- c("large", "medium", "small")
p <- c("low", "medium", "high")
type <- c("cauchy", "cauchyAff", "cauchy3d", "tdistDeg4", "tdistDeg2", "clover")

arg <- expand.grid(n = n, k = k, p = p, gamma = gamma, seed = seed)
arg <- sprintf("plot.R --n %d --k %d --p %s --gamma %d --seed %d",
                arg$n, arg$k, arg$p, arg$gamma, arg$seed)

connection <- file("arg.txt")
writeLines(arg, connection)
close(connection)
