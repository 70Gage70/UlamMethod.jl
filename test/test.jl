using UlamMethod
import Random
Random.seed!(1234)

### 1d
boundary1d = Boundary(0, 1)
traj1d = Trajectories(rand(1, 10), rand(1, 10))

ur = ulam_method(traj1d, LineBinner(10, boundary1d))

### 2d
n_points = 100
x0_rand = randn(2, n_points) + [fill(1, n_points) ;; fill(4, n_points)]'
xT_rand = x0_rand + rand(2, n_points)
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
n_points = 100000
x0_rand = randn(3, n_points) + [fill(1, n_points) ;; fill(1, n_points) ;; fill(4, n_points)]'
xT_rand = x0_rand + randn(3, n_points)
traj3d = Trajectories(x0_rand, xT_rand)

boundary3d = Boundary((0,0,0), (5, 5, 5))

ur = ulam_method(traj3d, HyperRectangleBinner(1000, boundary3d))
ur = ulam_method(traj3d, HyperRectangleBinner(1000, boundary3d), reinj_algo = SourceReinjection([(1, 2, 3)])) 
