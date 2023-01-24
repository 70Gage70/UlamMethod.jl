using UlamMethod
using Test
using HDF5

include("testdata.jl")

MATin = joinpath(@__DIR__, "x0x5-NA-undrogued.mat")
h5in = joinpath(@__DIR__, "x0x5-NA-undrogued.h5")

@testset "ulam_tpt_h5" begin
    for case in test_cases
        type = case[1]
        n_polys = case[2]
        f_ulam = joinpath(@__DIR__, "ulam_" * string(type) * "_" * string(n_polys) * ".h5")
        f_tpt = joinpath(@__DIR__, "TPT_" * string(type) * "_" * string(n_polys) * ".h5")
        @test ulam_method_tpt_test(h5in, n_polys, type, corners, A_centers, B_centers, f_ulam, f_tpt)
    end
end

@testset "ulam_tpt_MAT" begin
    for case in test_cases
        type = case[1]
        n_polys = case[2]
        f_ulam = joinpath(@__DIR__, "ulam_" * string(type) * "_" * string(n_polys) * ".h5")
        f_tpt = joinpath(@__DIR__, "TPT_" * string(type) * "_" * string(n_polys) * ".h5")
        @test ulam_method_tpt_test(MATin, n_polys, type, corners, A_centers, B_centers, f_ulam, f_tpt)
    end
end