module UlamMethod

using Meshes
using ArgCheck
using PolygonInbounds
using Graphs: SimpleDiGraph, strongly_connected_components
using ParallelKMeans
using PrecompileTools: @compile_workload 

include("earth-polygons.jl") # EarthPolygons module

include("hyperrectangle.jl")
export HyperRectangle

include("traj.jl")
export Trajectories

include("boundary.jl")
export Boundary, points

include("bins.jl")
export Bins, BinningAlgorithm

include(joinpath(@__DIR__, "..", "src", "binners", "membership-1d.jl"))
include(joinpath(@__DIR__, "..", "src", "binners", "membership-2d.jl"))

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
export HyperRectangleBinner

include("reinjection.jl")
export ReinjectionAlgorithm
export DataReinjection, SourceReinjection

include("result.jl")
export UlamResult
export P_open, P_closed, bins, bins_dis, membership

include("main.jl")
export ulam_method

@compile_workload begin
    import Random
    Random.seed!(1234)

    ### 1d
    boundary1d = Boundary(0, 1)
    traj1d = Trajectories(rand(1, 10), rand(1, 10))

    ur = ulam_method(traj1d, LineBinner(10, boundary1d))

    ### 2d
    n_points = 1000
    x0_rand = randn(2, n_points) + [fill(1, n_points) ;; fill(4, n_points)]'
    xT_rand = x0_rand + randn(2, n_points)
    traj2d = Trajectories(x0_rand, xT_rand)

    boundary2d = Boundary(0, 6, 0, 4)
    boundary2d = Boundary([(0,0),(6,0),(1,7),(1,6)])

    ur = ulam_method(traj2d, RectangleBinner(10, boundary2d))
    ur = ulam_method(traj2d, RectangleBinner(10, boundary2d, hardclip = false))
    ur = ulam_method(traj2d, RectangleBinner(10, boundary2d), reinj_algo = SourceReinjection([(1, 3)]))
    ur = ulam_method(traj2d, VoronoiBinner(10, boundary2d, traj2d))
    ur = ulam_method(traj2d, VoronoiBinner(10, boundary2d, traj2d), reinj_algo = SourceReinjection([(1, 3)]))
    ur = ulam_method(traj2d, TriangleBinner(10, boundary2d))
    ur = ulam_method(traj2d, TriangleBinner(10, boundary2d), reinj_algo = SourceReinjection([(1, 3)]))
    ur = ulam_method(traj2d, HexagonBinner(10, boundary2d))
    ur = ulam_method(traj2d, HexagonBinner(10, boundary2d), reinj_algo = SourceReinjection([(1, 3)]))

    ### ≥3D
    n_points = 1000
    x0_rand = randn(3, n_points) + [fill(1, n_points) ;; fill(1, n_points) ;; fill(4, n_points)]'
    xT_rand = x0_rand + randn(3, n_points)
    traj3d = Trajectories(x0_rand, xT_rand)
    
    boundary3d = Boundary((0,0,0), (5, 5, 5))
    
    ur = ulam_method(traj3d, HyperRectangleBinner(100, boundary3d))
    ur = ulam_method(traj3d, HyperRectangleBinner(100, boundary3d), reinj_algo = SourceReinjection([(1, 2, 3)]))
end

end # module