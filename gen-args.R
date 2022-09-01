# Make a txt file of simulation arguments, one argument combination per row
n <- 5000
k <- 400
p <- c("low", "medium", "high")
gamma <- 1
seed <- 204
path <- "figures/"

arg <- expand.grid(n = n, k = k, p = p, gamma = gamma, seed = seed,
                   path = path)
arg <- sprintf("plot.R --n %d --k %d --p %s --gamma %d --seed %d --path %s",
                arg$n, arg$k, arg$p, arg$gamma, arg$seed, arg$path)

connection <- file("arg.txt")
writeLines(arg, connection)
close(connection)
