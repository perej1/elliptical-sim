# Make a txt file of simulation arguments, one argument combination per row
n <- 5000
k <- 400
p <- c("low", "medium", "high")
gamma <- 1
seed <- 204

arg <- expand.grid(n = n, k = k, p = p, gamma = gamma, seed = seed)
arg <- sprintf("plot.R --n %d --k %d --p %s --gamma %d --seed %d",
                arg$n, arg$k, arg$p, arg$gamma, arg$seed)

connection <- file("arg.txt")
writeLines(arg, connection)
close(connection)
