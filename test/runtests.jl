using UlamMethod
using Test

@testset "UlamMethod.jl" begin
    @test UlamMethod.ulamTPT("x0x5-NA-undrogued", 500, "reg", extra_suffix = "_test")
end
