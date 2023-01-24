using UlamMethod
using Test
using HDF5


inpath = joinpath(@__DIR__, "x0x5-NA-undrogued")
regpath = joinpath(@__DIR__, "ulamTPT_reg_500.h5")
hexpath = joinpath(@__DIR__, "ulamTPT_hex_500.h5")
vorpath = joinpath(@__DIR__, "ulamTPT_vor_50.h5")

@testset "UlamMethod.jl" begin
    @test ulamTPTtest(inpath, 500, "reg", regpath) == true
    @test ulamTPTtest(inpath, 500, "hex", hexpath) == true
    @test ulamTPTtest(inpath, 50, "vor", vorpath) == true
end