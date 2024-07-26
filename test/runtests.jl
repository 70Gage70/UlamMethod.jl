using UlamMethod
using Test
import Random
Random.seed!(1234)

@testset "1D" begin
    traj1d = Trajectories(1, 1000)

    boundary1d = Boundary(0, 1)
    
    ur = ulam_method(traj1d, LineBinner(10, boundary1d))
end

@testset "2D" begin
    traj2d = Trajectories(2, 1000)

    boundary2d = Boundary([(0,0),(6,0),(1,7),(1,6)])
    boundary2d = Boundary(0, 6, 0, 4)

    ur = ulam_method(traj2d, RectangleBinner(10, boundary2d))
    ur = ulam_method(traj2d, RectangleBinner(10, boundary2d, hardclip = false))
    ur = ulam_method(traj2d, RectangleBinner(10, boundary2d), reinj_algo = SourceReinjection([(1, 3)]))
    ur = ulam_method(traj2d, VoronoiBinner(10, boundary2d, traj2d))
    ur = ulam_method(traj2d, VoronoiBinner(10, boundary2d, traj2d), reinj_algo = SourceReinjection([(1, 3)]))
    ur = ulam_method(traj2d, TriangleBinner(10, boundary2d))
    ur = ulam_method(traj2d, TriangleBinner(10, boundary2d), reinj_algo = SourceReinjection([(1, 3)]))
    ur = ulam_method(traj2d, HexagonBinner(10, boundary2d))
    ur = ulam_method(traj2d, HexagonBinner(10, boundary2d), reinj_algo = SourceReinjection([(1, 3)]))
    ur = ulam_method(traj2d, VoronoiBinner(100, boundary2d, traj2d), reinj_algo = SourceReinjection([(1, 3)]))
end

@testset "3D" begin
    traj3d = Trajectories(3, 100000)

    boundary3d = Boundary((0,0,0), (5, 5, 5))

    ur = ulam_method(traj3d, HyperRectangleBinner(1000, boundary3d))
    ur = ulam_method(traj3d, HyperRectangleBinner(1000, boundary3d), reinj_algo = SourceReinjection([(1, 2, 3)])) 
end
