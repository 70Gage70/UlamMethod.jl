using Meshes, CairoMakie
using ArgCheck
using PolygonInbounds
# coords(point).x to get values

include("bins.jl")
include("traj.jl")


### 1d
boundary1d = Segment((0,), (1,))
test_bin1d = bin(boundary1d, LineBinner(500))
traj1d = Trajectories(rand(1, 10), rand(1, 10))

### 2d
boundary2d = Ngon([(0,0),(6,0),(1,7),(1,6)]...)
test_bin2d = bin(boundary2d, RectangleBinner(500))
traj2d = Trajectories(rand(2, 10), rand(2, 10))

