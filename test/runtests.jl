using UlamMethod
using Test
using Meshes: sample, HomogeneousSampling, coords, Box
import Random
Random.seed!(1234)

@testset "1D" begin
    traj1d = Trajectories(1, 1000)

    boundary1d = AutoBoundary(traj1d)
    
    ur = ulam_method(traj1d, LineBinner(10, boundary1d))
end

@testset "2D Square" begin
    traj2d = Trajectories(2, 1000)

    boundary2d = Boundary(0, 6, 0, 4)
    pts = sample(boundary2d.boundary, HomogeneousSampling(100)) |> collect # for SourceReinjection
    pts = [(coords(pt).x.val, coords(pt).y.val) for pt in pts]

    ur = ulam_method(traj2d, RectangleBinner(100, boundary2d))
    ur = ulam_method(traj2d, RectangleBinner(100, boundary2d, hardclip = false))
    ur = ulam_method(traj2d, RectangleBinner(100, boundary2d), reinj_algo = SourceReinjection(pts))
    ur = ulam_method(traj2d, RectangleBinner(100, boundary2d), reinj_algo = StationaryReinjection())
    ur = ulam_method(traj2d, VoronoiBinner(10, boundary2d, traj2d))
    ur = ulam_method(traj2d, VoronoiBinner(10, boundary2d, traj2d), reinj_algo = SourceReinjection(pts))
    ur = ulam_method(traj2d, TriangleBinner(100, boundary2d))
    ur = ulam_method(traj2d, TriangleBinner(100, boundary2d), reinj_algo = SourceReinjection(pts))
    ur = ulam_method(traj2d, HexagonBinner(100, boundary2d))
    ur = ulam_method(traj2d, HexagonBinner(100, boundary2d), reinj_algo = SourceReinjection(pts))
    ur = ulam_method(traj2d, VoronoiBinner(10, boundary2d, traj2d), reinj_algo = SourceReinjection(pts))
end

@testset "2D Poly" begin
    traj2d = Trajectories(2, 1000)

    boundary2d = Boundary([(0,0),(6,0),(1,7),(1,6)])
    pts = sample(boundary2d.boundary, HomogeneousSampling(100)) |> collect # for SourceReinjection
    pts = [(coords(pt).x.val, coords(pt).y.val) for pt in pts]

    ur = ulam_method(traj2d, RectangleBinner(100, boundary2d))
    ur = ulam_method(traj2d, RectangleBinner(100, boundary2d, hardclip = false))
    ur = ulam_method(traj2d, RectangleBinner(100, boundary2d), reinj_algo = SourceReinjection(pts))
    ur = ulam_method(traj2d, RectangleBinner(100, boundary2d), reinj_algo = StationaryReinjection())
    ur = ulam_method(traj2d, VoronoiBinner(10, boundary2d, traj2d))
    ur = ulam_method(traj2d, VoronoiBinner(10, boundary2d, traj2d), reinj_algo = SourceReinjection(pts))
    ur = ulam_method(traj2d, TriangleBinner(100, boundary2d))
    ur = ulam_method(traj2d, TriangleBinner(100, boundary2d), reinj_algo = SourceReinjection(pts))
    ur = ulam_method(traj2d, HexagonBinner(100, boundary2d))
    ur = ulam_method(traj2d, HexagonBinner(100, boundary2d), reinj_algo = SourceReinjection(pts))
    ur = ulam_method(traj2d, VoronoiBinner(10, boundary2d, traj2d), reinj_algo = SourceReinjection(pts))
end

@testset "2D Auto Rec" begin
    traj2d = Trajectories(2, 1000)

    boundary2d = AutoBoundary(traj2d)
    pts = sample(boundary2d.boundary, HomogeneousSampling(100)) |> collect # for SourceReinjection
    pts = [(coords(pt).x.val, coords(pt).y.val) for pt in pts]

    ur = ulam_method(traj2d, RectangleBinner(100, boundary2d))
    ur = ulam_method(traj2d, RectangleBinner(100, boundary2d, hardclip = false))
    ur = ulam_method(traj2d, RectangleBinner(100, boundary2d), reinj_algo = SourceReinjection(pts))
    ur = ulam_method(traj2d, RectangleBinner(100, boundary2d), reinj_algo = StationaryReinjection())
    ur = ulam_method(traj2d, VoronoiBinner(10, boundary2d, traj2d))
    ur = ulam_method(traj2d, VoronoiBinner(10, boundary2d, traj2d), reinj_algo = SourceReinjection(pts))
    ur = ulam_method(traj2d, TriangleBinner(100, boundary2d))
    ur = ulam_method(traj2d, TriangleBinner(100, boundary2d), reinj_algo = SourceReinjection(pts))
    ur = ulam_method(traj2d, HexagonBinner(100, boundary2d))
    ur = ulam_method(traj2d, HexagonBinner(100, boundary2d), reinj_algo = SourceReinjection(pts))
    ur = ulam_method(traj2d, VoronoiBinner(100, boundary2d, traj2d), reinj_algo = SourceReinjection(pts))
end

@testset "2D Auto Poly" begin
    traj2d = Trajectories(2, 1000)

    boundary2d = AutoBoundary2D(traj2d)
    pts = sample(boundary2d.boundary, HomogeneousSampling(100)) |> collect # for SourceReinjection
    pts = [(coords(pt).x.val, coords(pt).y.val) for pt in pts]

    ur = ulam_method(traj2d, RectangleBinner(100, boundary2d))
    ur = ulam_method(traj2d, RectangleBinner(100, boundary2d, hardclip = false))
    ur = ulam_method(traj2d, RectangleBinner(100, boundary2d), reinj_algo = SourceReinjection(pts))
    ur = ulam_method(traj2d, RectangleBinner(100, boundary2d), reinj_algo = StationaryReinjection())
    ur = ulam_method(traj2d, VoronoiBinner(10, boundary2d, traj2d))
    ur = ulam_method(traj2d, VoronoiBinner(10, boundary2d, traj2d), reinj_algo = SourceReinjection(pts))
    ur = ulam_method(traj2d, TriangleBinner(100, boundary2d))
    ur = ulam_method(traj2d, TriangleBinner(100, boundary2d), reinj_algo = SourceReinjection(pts))
    ur = ulam_method(traj2d, HexagonBinner(100, boundary2d))
    ur = ulam_method(traj2d, HexagonBinner(100, boundary2d), reinj_algo = SourceReinjection(pts))
    ur = ulam_method(traj2d, VoronoiBinner(100, boundary2d, traj2d), reinj_algo = SourceReinjection(pts))
end

@testset "3D Box" begin
    traj3d = Trajectories(3, 100000)

    boundary3d = Boundary((0,0,0), (5, 5, 5))
    box = Box(boundary3d.boundary.min, boundary3d.boundary.max)
    pts = sample(box, HomogeneousSampling(100)) |> collect # for SourceReinjection
    pts = [(coords(pt).x.val, coords(pt).y.val, coords(pt).z.val) for pt in pts]


    ur = ulam_method(traj3d, HyperRectangleBinner(1000, boundary3d))
    ur = ulam_method(traj3d, HyperRectangleBinner(1000, boundary3d), reinj_algo = SourceReinjection(pts))
    ur = ulam_method(traj3d, HyperRectangleBinner(1000, boundary3d), reinj_algo = StationaryReinjection()) 
end

@testset "3D Auto" begin
    traj3d = Trajectories(3, 100000)

    boundary3d = AutoBoundary(traj3d)
    box = Box(boundary3d.boundary.min, boundary3d.boundary.max)
    pts = sample(box, HomogeneousSampling(100)) |> collect # for SourceReinjection
    pts = [(coords(pt).x.val, coords(pt).y.val, coords(pt).z.val) for pt in pts]

    ur = ulam_method(traj3d, HyperRectangleBinner(1000, boundary3d))
    ur = ulam_method(traj3d, HyperRectangleBinner(1000, boundary3d), reinj_algo = SourceReinjection(pts))
    ur = ulam_method(traj3d, HyperRectangleBinner(1000, boundary3d), reinj_algo = StationaryReinjection())  
end

@testset "4D" begin
    traj4d = Trajectories(4, 100000)

    boundary4d = AutoBoundary(traj4d)

    ur = ulam_method(traj4d, HyperRectangleBinner(1000, boundary4d))
end