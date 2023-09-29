# Generate arguments for simulation scenarios

type <- c("cauchy", "cauchyAff", "tdistDeg4", "tdistSkew")
s <- c(1, 100)
d <- 2:3
n <- c(1000, 5000)
p <- c("low", "medium", "high")
k <- c("small", "medium", "large")
seed <- 278

arg_2d <- expand.grid(
  type = type[-4],
  s = s[2],
  d = d[1],
  n = n,
  p = p,
  k = k,
  seed = seed
)

arg_3d <- expand.grid(
  type = type[c(1, 3)],
  s = s[2],
  d = d[2],
  n = n,
  p = p,
  k = k,
  seed = seed
)

arg_skew <- expand.grid(
  type = type[4],
  s = s[1],
  d = d[1],
  n = n,
  p = p,
  k = k,
  seed = seed
)

arg <- rbind(arg_3d, arg_2d, arg_skew)
arg_vector <- sprintf(paste0("simulate.R --type %s --s %d --d %d --n %d ",
                      "--p %s --k %s --seed %d"),
  arg$type,
  arg$s,
  arg$d,
  arg$n,
  arg$p,
  arg$k,
  arg$seed
)

readr::write_lines(arg_vector, "sim-args.txt")
readr::write_csv(arg, "sim-args.csv")
