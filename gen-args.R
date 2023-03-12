# Simulation arguments
#type <- c("cauchy", "cauchyAff")
type <- c("cauchy", "cauchyAff", "tdistDeg4", "tdistDeg2", "clover")
n <- c(500, 1000, 5000)
seed <- 204
p <- c("low", "medium", "high")
k <-c("large", "medium", "small")
i <- 1

# Arguments for creating samples
arg <- expand.grid(type = type, n = n, seed = seed)
arg <- sprintf("gen-sample.R --type %s --n %d --seed %d",
               arg$type,
               arg$n,
               arg$seed)

connection <- file("args/arg-gen-sample.txt")
writeLines(arg, connection)
close(connection)

# # Arguments for estimation
# arg <- expand.grid(type = type, n = n, p = p, k = k)
# arg <- sprintf("estimate.R --type %s --n %d --p %s --k %s",
#                arg$type,
#                arg$n,
#                arg$p,
#                arg$k)
# 
# connection <- file("args/arg-estimate.txt")
# writeLines(arg, connection)
# close(connection)
# 
# # Arguments for plotting
# arg <- expand.grid(type = type, n = n, p = p, k = k, i = 1)
# arg <- sprintf("plot.R --type %s --n %d --p %s --k %s --i %d",
#                arg$type,
#                arg$n,
#                arg$p,
#                arg$k,
#                arg$i)
# 
# connection <- file("args/arg-plot.txt")
# writeLines(arg, connection)
# close(connection)
