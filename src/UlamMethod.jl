module UlamMethod

using Meshes
using ArgCheck
using PolygonInbounds
using Graphs: SimpleDiGraph, strongly_connected_components
using ParallelKMeans
using PrecompileTools: @compile_workload 

include("earth-polygons.jl") # EarthPolygons module

include("traj.jl")
export Trajectories

include("boundary.jl")
export Boundary, points

include("bins.jl")
export Bins, BinningAlgorithm
export LineBinner, RectangleBinner, TriangleBinner, HexagonBinner, VoronoiBinner
export bin

include("membership.jl")
export membership

include("reinjection.jl")
export ReinjectionAlgorithm
export DataReinjection, SourceReinjection
export reinject!

include("result.jl")
export UlamResult
export P_open, P_closed, bins, bins_dis

include("main.jl")
export ulam_method

@compile_workload begin
    import Random
    Random.seed!(1234)

    ### 1d
    boundary1d = Boundary(0, 1)
    traj1d = Trajectories(rand(1, 10), rand(1, 10))

    ur = ulam_method(traj1d, boundary1d, LineBinner(10))

    ### 2d
    n_points = 100
    x0_rand = randn(2, n_points) + [fill(1, n_points) ;; fill(4, n_points)]'
    xT_rand = x0_rand + rand(2, n_points)
    traj2d = Trajectories(x0_rand, xT_rand)

    boundary2d = Boundary(0, 6, 0, 4)
    boundary2d = Boundary([(0,0),(6,0),(1,7),(1,6)])

    ur = ulam_method(traj2d, boundary2d, RectangleBinner(10))
    ur = ulam_method(traj2d, boundary2d, RectangleBinner(10), reinj_algo = SourceReinjection([(1, 3)]))
    ur = ulam_method(traj2d, boundary2d, VoronoiBinner(10, traj2d))
    ur = ulam_method(traj2d, boundary2d, VoronoiBinner(10, traj2d), reinj_algo = SourceReinjection([(1, 3)]))
    ur = ulam_method(traj2d, boundary2d, TriangleBinner(10))
    ur = ulam_method(traj2d, boundary2d, TriangleBinner(10), reinj_algo = SourceReinjection([(1, 3)]))
    ur = ulam_method(traj2d, boundary2d, HexagonBinner(10))
    ur = ulam_method(traj2d, boundary2d, HexagonBinner(10), reinj_algo = SourceReinjection([(1, 3)]))
end

end # module