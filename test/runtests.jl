using Test
using EnergyStatistics

@testset "Double Center" begin
    x = [1, 1, -1, -1]
    dm = Array{Float64}([0 0 2 2; 0 0 2 2; 2 2 0 0; 2 2 0 0])
    @test EnergyStatistics.DistanceMatrix(x) == dm

    dcdm = Array{Float64}([-1 -1 1 1; -1 -1 1 1; 1 1 -1 -1; 1 1 -1 -1])
    @test EnergyStatistics.dcenter!(dm) == dcdm
    @test EnergyStatistics.dcenter!(EnergyStatistics.DistanceMatrix(x)) == dcdm

    x = [-1, 0, 1]
    A = EnergyStatistics.DistanceMatrix(x)
    dm = Array{Float64}([0 1 2; 1 0 1; 2 1 0])
    @test A == dm

    EnergyStatistics.ucenter!(A)
    @test A == Array{Float64}([0 0 0; 0 0 0; 0 0 0])


end

@testset "Matrix Algebra" begin
    L = 1000
    x = randn(L)
    y = 0.1 * randn(L) .+ x

    A = EnergyStatistics.dcenter!(EnergyStatistics.DistanceMatrix(x))
    B = EnergyStatistics.dcenter!(EnergyStatistics.DistanceMatrix(y))

    @test dcov(A, B) == dcov(x,y)
    @test dvar(A)    == dvar(x)
    @test dcor(A, B) == dcor(x,y)

    dvarA = dvar(A)
    dvarB = dvar(B)
    @test dcor(A, B) == dcor(A, B, dvarA, dvarB)
end

@testset "Values" begin
    x = [1.0,  1.0, -1.0, -1.0]
    y = [1.0, -1.0,  1.0, -1.0]
    @test dcor(x,y) == 0.0

    x = [1,  1, -1, -1]
    y = [1, -1,  1, -1]
    @test dcor(x,y) == 0.0

    x = [-1.0, -1.0,  1.0,  1.0]
    @test dcor(x,x) == 1.0

    x = [-1.0, -1.0,  1.0,  1.0]
    y = [1.0,   1.0,  1.0,  1.0]
    @test dcov(x,y) == 0.0

    x = [-1.0, 0.0,  1.0]
    y = [1.0,  0.0,  1.0]
    @test dcor(x,y) ≈ 0.5623413251903491

    x = collect(-1:0.01:1)
    y = @. x^4 - x^2
    @test dcor(x, y) ≈ 0.3742040504583154

end

;
