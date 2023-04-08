using .UlamTypes
using BenchmarkTools

include("main.jl")

fin = "x0x5-NA-undrogued.mat"
corners = [-100, 15, -9, 39]
poly_number = 500
poly_type = "hex"
traj = UlamTrajectories(fin);
# domain = UlamDomain(corners..., poly_number = poly_number, poly_type = poly_type);
domain = UlamDomain(corners..., poly_type = "vor")
ulam = ulam_method(traj, domain)


function utest()
    traj = UlamTrajectories(fin);
    domain = UlamDomain(corners..., poly_number = poly_number, poly_type = poly_type);
    ulam = ulam_method(traj, domain)

    return ulam
end