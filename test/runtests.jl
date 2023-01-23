using UlamMethod
# include("../src/UlamMethod.jl")
using Test

@testset "UlamMethod.jl" begin
    @test UlamMethod.ulamTPTtest("x0x5-NA-undrogued", 500, "reg", "ulamTPT_reg_500.h5") == true
    @test UlamMethod.ulamTPTtest("x0x5-NA-undrogued", 500, "hex", "ulamTPT_hex_500.h5") == true
    @test UlamMethod.ulamTPTtest("x0x5-NA-undrogued", 50, "vor", "ulamTPT_vor_50.h5") == true
end