using .UlamTypes
using BenchmarkTools

include("main.jl")

fin = "data/x0x5-NA-undrogued.mat"
corners = [-100, 15, -9, 39]
poly_number = 50
poly_type = "vor"
traj = UlamTrajectories(fin);
domain = UlamDomain(corners..., poly_number = poly_number, poly_type = poly_type);
ulam = ulam_method(traj, domain)


function utest()
    traj = UlamTrajectories(fin);
    domain = UlamDomain(corners..., poly_number = poly_number, poly_type = poly_type);
    ulam = ulam_method(traj, domain)

    return ulam
end