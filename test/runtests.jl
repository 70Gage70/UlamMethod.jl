using Test
using UlamMethod, HDF5

include("test-functions.jl")
include("test-cases.jl")

@testset "ulam_method" begin
    for t in test_cases
        traj, domain, ftest = t
        res = ulam_method(traj, domain)
        @test ulam_test(ftest, res)

    end
end
