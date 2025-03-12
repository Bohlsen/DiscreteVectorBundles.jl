using DiscreteVectorBundles
using Test

@testset "DiscreteVectorBundles.jl" begin

    include("testeigenbundles.jl")

    @test mobiusbundletest() == 1
    @test spherebundletest() == 2

    @test geofluidstest(1) == 1
    @test geofluidstest(2) == 0
    @test geofluidstest(3) == 2
end
