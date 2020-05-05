# These benchmarks are not included in the runtests.jl to avoid
# heavy dependencies on RCall, R and the energy package in R.

using BenchmarkTools
using EnergyStatistics
using RCall


N = 4272
x = rand(N)
y = rand(N)

@benchmark dcor($x, $y)


time_in_milliseconds = convert(Float64, R"
library(energy)
x <- runif($N)
y <- runif($N)

# not sure whether this is the best way to benchmark R code
start_time <- Sys.time()
dcor(x, y)
end_time <- Sys.time()

as.numeric(end_time - start_time) * 1000
")

# on my MacBook julia is about 20x faster than R, which is supposed to
# utilize a c implementation.
