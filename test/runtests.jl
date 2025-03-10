using DiscreteVectorBundles
using Test

@testset "DiscreteVectorBundles.jl" begin

    include("testeigenbundles.jl")
    @test spherebundletest() == 2
end
