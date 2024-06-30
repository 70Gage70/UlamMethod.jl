using UlamMethod
using Distributions # FOR TESTING ONLY


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

ur_rec = ulam_method(traj2d, boundary2d, RectangleBinner(500), reinj_algo = SourceReinjection([(1, 3)]))
ur_vor = ulam_method(traj2d, boundary2d, VoronoiBinner(500, traj2d), reinj_algo = SourceReinjection([(1, 3)]))
ur_hex = ulam_method(traj2d, boundary2d, HexagonBinner(500), reinj_algo = SourceReinjection([(1, 3)]))
ur_tri = ulam_method(traj2d, boundary2d, TriangleBinner(500), reinj_algo = SourceReinjection([(1, 3)]))
