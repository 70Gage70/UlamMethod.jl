module UlamMethod

using Meshes
using ArgCheck
using PolygonInbounds
using Graphs: SimpleDiGraph, strongly_connected_components
using ParallelKMeans
import Distributions
using StatsBase: countmap
using PrecompileTools: @compile_workload 

include("traj.jl")
export Trajectories

include("bins.jl")
export Bins, BinningAlgorithm

include(joinpath(@__DIR__, "..", "src", "binners", "membership-1d.jl"))
include(joinpath(@__DIR__, "..", "src", "binners", "membership-2d.jl"))

include("boundary.jl")
export Boundary, points, AutoBoundary, AutoBoundary2D

include(joinpath(@__DIR__, "..", "src", "binners", "line.jl"))
export LineBinner

include(joinpath(@__DIR__, "..", "src", "binners", "rectangle.jl"))
export RectangleBinner

include(joinpath(@__DIR__, "..", "src", "binners", "triangle.jl"))
export TriangleBinner

include(joinpath(@__DIR__, "..", "src", "binners", "hexagon.jl"))
export HexagonBinner

include(joinpath(@__DIR__, "..", "src", "binners", "voronoi.jl"))
export VoronoiBinner

include(joinpath(@__DIR__, "..", "src", "binners", "hyperrectangle.jl"))
export HyperRectangle, HyperRectangleBinner

include("reinjection.jl")
export ReinjectionAlgorithm
export DataReinjection, SourceReinjection

include("result.jl")
export UlamResult
export P_open, P_closed, bins, bins_dis, membership

include("main.jl")
export ulam_method

include("show.jl")
export show

include("viz.jl")
export viz, viz!

include("earth-polygons.jl") # EarthPolygons module

@compile_workload begin
    import Random
    Random.seed!(1234)

    ### 1d
    traj1d = Trajectories(1, 1000)

    boundary2d = AutoBoundary(traj1d)
    boundary1d = Boundary(0, 1)
    
    ur = ulam_method(traj1d, LineBinner(100, boundary1d))

    ### 2d
    traj2d = Trajectories(2, 1000)

    boundary2d = AutoBoundary(traj2d)
    boundary2d = AutoBoundary2D(traj2d)
    boundary2d = Boundary([(0,0),(6,0),(1,7),(1,6)])
    boundary2d = Boundary(0, 6, 0, 4)

    ur = ulam_method(traj2d, RectangleBinner(100, boundary2d))
    ur = ulam_method(traj2d, RectangleBinner(100, boundary2d, hardclip = false))
    ur = ulam_method(traj2d, RectangleBinner(100, boundary2d), reinj_algo = SourceReinjection([(1, 3)]))
    ur = ulam_method(traj2d, VoronoiBinner(10, boundary2d, traj2d))
    ur = ulam_method(traj2d, VoronoiBinner(10, boundary2d, traj2d), reinj_algo = SourceReinjection([(1, 3)]))
    ur = ulam_method(traj2d, TriangleBinner(100, boundary2d))
    ur = ulam_method(traj2d, TriangleBinner(100, boundary2d), reinj_algo = SourceReinjection([(1, 3)]))
    ur = ulam_method(traj2d, HexagonBinner(100, boundary2d))
    ur = ulam_method(traj2d, HexagonBinner(100, boundary2d), reinj_algo = SourceReinjection([(1, 3)]))
    ur = ulam_method(traj2d, VoronoiBinner(10, boundary2d, traj2d), reinj_algo = SourceReinjection([(1, 3)]))

    ### ≥3D
    traj3d = Trajectories(3, 100000)

    boundary3d = AutoBoundary(traj3d)
    boundary3d = Boundary((0,0,0), (5, 5, 5))

    ur = ulam_method(traj3d, HyperRectangleBinner(1000, boundary3d))
    ur = ulam_method(traj3d, HyperRectangleBinner(1000, boundary3d), reinj_algo = SourceReinjection([(1, 2, 3)])) 
end

end # module