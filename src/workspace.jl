using .UlamTypes

include("ulam-binners.jl")
include("helpers.jl")

traj = UlamTrajectories("x0x5-NA-undrogued.mat");
domain = UlamDomain(-100, 15, -9, 39, poly_number = 500);
polys = square_binner(traj, domain);
