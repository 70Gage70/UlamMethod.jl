using Meshes
import CairoMakie as Mke
using ArgCheck
using PolygonInbounds
using Graphs: SimpleDiGraph, strongly_connected_components
using ParallelKMeans

using Distributions # FOR TESTING ONLY

include("traj.jl")
include("boundary.jl")
include("bins.jl")
include("membership.jl")
include("reinjection.jl")
include("result.jl")
include("main.jl")


### 1d
boundary1d = Boundary(0, 1)
traj1d = Trajectories(rand(1, 10), rand(1, 10))

ur = ulam_method(traj1d, boundary1d, LineBinner(500))

### 2dx
x0_rand = rand(MvNormal([1, 4], [1 0; 0 1]), 10000)
xT_rand = x0_rand + rand(MvNormal([0, 0], [1 0; 0 1]), size(x0_rand, 2))
traj2d = Trajectories(x0_rand, xT_rand)

boundary2d = Boundary([(0,0),(6,0),(1,7),(1,6)])
# boundary2d = Boundary(0, 6, 0, 4)

ur = ulam_method(traj2d, boundary2d, RectangleBinner(500), reinj_algo = SourceReinjection([(1, 3)]))
ur = ulam_method(traj2d, boundary2d, VoronoiBinner(500, traj2d), reinj_algo = SourceReinjection([(1, 3)]))