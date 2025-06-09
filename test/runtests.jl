using DiscreteVectorBundles
using Test

@testset "DiscreteVectorBundles.jl" begin

    include("testeigenbundles.jl")

    @test mobiusbundletest() == 1
    @test spherebundletest() == 2

    #Signs here assigned to agree with Delplace's lecture notes
    #where the expected spectral flow index is -1*C1 (Delplace does somewhat explain this)
    @test diractest()[1] == 1
    @test diractest()[2] == -1

    #Signs here (and the operator H_geophys) chosen to agree with Delplace's lecture notes
    #it does seem to depend sensitively on the choices of the definition of the spectral flow
    @test geofluidstest(1) == 2
    @test geofluidstest(2) == 0
    @test geofluidstest(3) == -2

    @test TLCWtest(1) == 0
    @test TLCWtest(2) == 0
    @test TLCWtest(3) == -1
    @test TLCWtest(4) == 1
    @test TLCWtest(5) == 0
    @test TLCWtest(6) == -1
    @test TLCWtest(7) == 1
    @test TLCWtest(8) == 0
    @test TLCWtest(9) == 0

    @test dirac2DWignermatrixsampletest() == -1
end
