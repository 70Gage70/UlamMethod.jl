module UlamMethod

using Meshes
using ArgCheck
using PolygonInbounds
using Graphs: SimpleDiGraph, strongly_connected_components
using ParallelKMeans
import Distributions
using StatsBase: countmap
using LinearAlgebra: eigvecs
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
export DataReinjection, UniformReinjection, SourceReinjection, StationaryReinjection

include("result.jl")
export UlamResult
export P_open, P_closed, bins, bins_dis, membership, counts

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
    boundary2d = Boundary([(0,0),(6,0),(1,7),(1,6)])
    boundary2d = Boundary(0, 6, 0, 4)

    boundary2d = AutoBoundary2D(traj2d)
    pts = sample(boundary2d.boundary, HomogeneousSampling(100)) |> collect # for SourceReinjection
    pts = [(coords(pt).x.val, coords(pt).y.val) for pt in pts]

    ur = ulam_method(traj2d, RectangleBinner(100, boundary2d))
    ur = ulam_method(traj2d, RectangleBinner(100, boundary2d, hardclip = false))
    ur = ulam_method(traj2d, RectangleBinner(100, boundary2d), reinj_algo = UniformReinjection())
    ur = ulam_method(traj2d, RectangleBinner(100, boundary2d), reinj_algo = SourceReinjection(pts))
    ur = ulam_method(traj2d, RectangleBinner(100, boundary2d), reinj_algo = StationaryReinjection())
    ur = ulam_method(traj2d, VoronoiBinner(10, boundary2d, traj2d))
    ur = ulam_method(traj2d, VoronoiBinner(10, boundary2d, traj2d), reinj_algo = SourceReinjection(pts))
    ur = ulam_method(traj2d, TriangleBinner(100, boundary2d))
    ur = ulam_method(traj2d, TriangleBinner(100, boundary2d), reinj_algo = SourceReinjection(pts))
    ur = ulam_method(traj2d, HexagonBinner(100, boundary2d))
    ur = ulam_method(traj2d, HexagonBinner(100, boundary2d), reinj_algo = SourceReinjection(pts))
    ur = ulam_method(traj2d, VoronoiBinner(100, boundary2d, traj2d), reinj_algo = SourceReinjection(pts))

    ### ≥3D
    traj3d = Trajectories(3, 100000)

    boundary3d = AutoBoundary(traj3d)
    boundary3d = Boundary((0,0,0), (5, 5, 5))

    ur = ulam_method(traj3d, HyperRectangleBinner(1000, boundary3d))
    ur = ulam_method(traj3d, HyperRectangleBinner(1000, boundary3d), reinj_algo = SourceReinjection([(1, 2, 3)])) 
end

using MAT

function _ulam(ARGS_bin_type::String)
  ######### EDITING SECTION #########
  ###################################

  # path to the mat file with data and variable names
  file_in = "/Users/gagebonner/Desktop/Repositories/SparseSargassum.jl/data/x0x5-2018-water.mat"
  x0_name, xT_name, y0_name, yT_name = "x0", "xT", "y0", "yT"

  bin_type = ARGS_bin_type    # bin type, should be one of "rec", "tri", "hex" or "vor"
  n_bins = 500        # number of bins requested

  # The fraction of data to place in nirvana on each side of the computational domain
  # e.g. if ω_bot = 0.02, then the box will be drawn such that 2% of the data are below the box
  ω_left  = 0.00
  ω_right = 0.00
  ω_bot   = 0.02
  ω_top   = 0.02

  # reinjection algorithm, should be one of "data", "stat" or "unif"
  # "data": reinject according to actual trajectories
  # "stat": reinject according to the stationary distribution
  # "unif": reinject uniformly
  reinj_type = "data"

  # output file
  # the output will be created here with the following 
  # the output contains three fields, "P", "bins_verts" and "info"
  # - P is P_closed
  # - bins_verts is a list of vertices such that the ith entry in bins_verts corresponds to 
  #   the ith row/column of P. The vertices are N_points x 2 matrices, e.g. they will be 
  #   4 x 2 matrices for rectangular bins.
  # - info is this information
  file_out = "/Users/gagebonner/Desktop/Repositories/SparseSargassum.jl/data/ulam-out.mat"




  ########## EDIT PAST HERE AT YOUR OWN RISK :) ##########
  ########################################################
  ########################################################

  x0, xT, y0, yT = matread(file_in) |> f -> (f[x0_name], f[xT_name], f[y0_name], f[yT_name]) .|> vec
  traj = Trajectories(permutedims(hcat(x0, y0)), permutedims(hcat(xT, yT)))
  boundary = AutoBoundary(traj, nirvana = [(ω_left, ω_right), (ω_bot, ω_top)])

  if bin_type == "rec"
      binner = RectangleBinner(n_bins, boundary, hardclip = false)
  elseif bin_type == "tri"
      binner = TriangleBinner(n_bins, boundary, hardclip = false)
  elseif bin_type == "hex"
      binner = HexagonBinner(n_bins, boundary, hardclip = false)
  elseif bin_type == "vor"
      binner = VoronoiBinner(n_bins, boundary, traj, hardclip = false)
  else
      error("Unsupported bin type.")
  end

  if reinj_type == "data"
      reinj_algo = DataReinjection()
  elseif reinj_type == "stat"  
      reinj_algo = StationaryReinjection()
  elseif reinj_type == "unif"  
      reinj_algo = UniformReinjection()
  end

  ulam = ulam_method(traj, binner, reinj_algo = reinj_algo)
  P = P_closed(ulam)
  bins_verts = points(bins(ulam))
  bins_verts = [b .|> collect |> x -> stack(x, dims = 1) for b in bins_verts]
  info = "P is P_closed. bins_verts is a list of vertices such that the ith entry in bins_verts \
  corresponds to the ith row/column of P. The vertices are N_points x 2 matrices, e.g. they will be \
  4 x 2 matrices for rectangular bins."
  matwrite(file_out, Dict("P" => P, "bins_verts" => bins_verts, "info" => info))

  return nothing
end

function julia_main()::Cint
    _ulam(ARGS[1])
    return 0 # if things finished successfully
  end

end # module