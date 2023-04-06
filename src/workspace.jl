using .UlamTypes

include("ulam-binners.jl")

Ut = UlamTrajectories("x0x5-NA-undrogued.mat");
Ud = UlamDomain(-100, 15, -9, 39, bin_number = 500);
polys = square_binner(Ut, Ud);
