# These tests are not included in the runtests.jl to avoid
# heavy dependencies on RCall, R and the energy package in R.

using EnergyStatistics
using RCall
using Test

@testset "compare to R" begin
    N = 1000
    x = rand(N)
    y = rand(N);

    @test convert(Float64, R"library(energy); dcov($x, $y)") ≈ dcov(x,y)
    @test convert(Float64, R"library(energy); dcor($x, $y)") ≈ dcor(x,y)
    @test convert(Float64, R"library(energy); dcov($x, $x)") ≈ dvar(x)

    @test convert(Float64, R"library(energy); DCOR($x, $y)$dCov") ≈ dcov(x,y)
    @test convert(Float64, R"library(energy); DCOR($x, $y)$dCor") ≈ dcor(x,y)
    @test convert(Float64, R"library(energy); DCOR($x, $y)$dVarX") ≈ dvar(x)
    @test convert(Float64, R"library(energy); DCOR($x, $y)$dVarY") ≈ dvar(y)

    dx = EnergyStatistics.dcenter!(EnergyStatistics.DistanceMatrix(x));
    @test convert(Matrix{Float64}, R"library(energy); Dcenter(dist($x))") ≈ dx

    dx = EnergyStatistics.ucenter!(EnergyStatistics.DistanceMatrix(x));
    @test convert(Matrix{Float64}, R"library(energy); Ucenter(dist($x))") ≈ dx
end
