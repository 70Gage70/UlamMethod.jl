using .UlamTypes
using BenchmarkTools

include("main.jl")

function utest()
    traj = UlamTrajectories("x0x5-NA-undrogued.mat");
    domain = UlamDomain(-100, 15, -9, 39, poly_number = 100, poly_type = "hex");
    ulam = ulam_method(traj, domain)

    return ulam
end