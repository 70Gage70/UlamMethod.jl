using .UlamTypes
using Graphs:SimpleDiGraph,  strongly_connected_components

include("helpers.jl")
include("ulam-binner.jl")
include("ulam-nirvana.jl")


traj = UlamTrajectories("x0x5-NA-undrogued.mat");
domain = UlamDomain(-100, 15, -9, 39, poly_number = 500);
polys = binner_square(domain);

#########################################################################################################
#########################################################################################################
#########################################################################################################

ulam = ulam_nirvana(traj, domain, polys)