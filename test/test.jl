using UlamMethod
using Distributions # FOR TESTING ONLY


### 1d
boundary1d = Boundary(0, 1)
traj1d = Trajectories(rand(1, 10), rand(1, 10))

ur = ulam_method(traj1d, LineBinner(10, boundary1d))

### 2dx
x0_rand = rand(MvNormal([1, 4], [1 0; 0 1]), 10000)
xT_rand = x0_rand + rand(MvNormal([0, 0], [1 0; 0 1]), size(x0_rand, 2))
traj2d = Trajectories(x0_rand, xT_rand)

boundary2d = Boundary([(0,0),(6,0),(1,7),(1,6)])
# boundary2d = Boundary(0, 6, 0, 4)

### 2d
n_points = 100000
x0_rand = randn(2, n_points) + [fill(1, n_points) ;; fill(4, n_points)]'
xT_rand = x0_rand + rand(2, n_points)
traj2d = Trajectories(x0_rand, xT_rand)

boundary2d = Boundary(0, 6, 0, 4)
boundary2d = Boundary([(0,0),(6,0),(1,7),(1,6)])

ur_rec = ulam_method(traj2d, RectangleBinner(100, boundary2d))
ur = ulam_method(traj2d, RectangleBinner(100, boundary2d, hardclip = false))
ur = ulam_method(traj2d, RectangleBinner(100, boundary2d), reinj_algo = SourceReinjection([(1, 3)]))
ur = ulam_method(traj2d, VoronoiBinner(100, boundary2d, traj2d))
ur = ulam_method(traj2d, VoronoiBinner(100, boundary2d, traj2d), reinj_algo = SourceReinjection([(1, 3)]))
ur = ulam_method(traj2d, TriangleBinner(100, boundary2d))
ur = ulam_method(traj2d, TriangleBinner(100, boundary2d), reinj_algo = SourceReinjection([(1, 3)]))
ur = ulam_method(traj2d, HexagonBinner(100, boundary2d))
ur = ulam_method(traj2d, HexagonBinner(100, boundary2d), reinj_algo = SourceReinjection([(1, 3)]))
