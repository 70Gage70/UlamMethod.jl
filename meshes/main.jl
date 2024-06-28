using Meshes, CairoMakie
using ArgCheck
using PolygonInbounds

include("boundary.jl")
include("bins.jl")
include("traj.jl")
include("membership.jl")

### 1d
boundary1d = Segment((0,), (1,)) |> Boundary
bin1d = bin(boundary1d, LineBinner(500))
traj1d = Trajectories(rand(1, 10), rand(1, 10))
membs1d = membership(traj1d, bin1d)

### 2d
boundary2d = Ngon([(0,0),(6,0),(1,7),(1,6)]...) |> Boundary
bin2d = bin(boundary2d, RectangleBinner(500))
traj2d = Trajectories(5*(rand(2, 10) .- 0.5), 5*rand(2, 10))
membs2d = membership(traj2d, bin2d)

