using Meshes, CairoMakie
using ArgCheck
using PolygonInbounds

include("boundary.jl")
include("bins.jl")
include("traj.jl")

### 1d
boundary1d = Segment((0,), (1,)) |> Boundary
bin1d = bin(boundary1d, LineBinner(500))
traj1d = Trajectories(rand(1, 10), rand(1, 10))

### 2d
boundary2d = Ngon([(0,0),(6,0),(1,7),(1,6)]...) |> Boundary
bin2d = bin(boundary2d, RectangleBinner(500))
traj2d = Trajectories(rand(2, 10), rand(2, 10))

nodes = []
edges = []

for i = 1:length(bin2d.bins)
    bin_i = bin2d.bins[i]
    push!(nodes, collect(coords.(bin_i.vertices) .|> x -> [x.x.val, x.y.val]))
    node_idx = size(nodes[i], 1)
    push!(edges, [[k, k == node_idx ? 1 : k + 1, i] for k = 1:node_idx])
end

nodes = vcat(stack.(nodes, dims = 1)...)
edges = vcat(stack.(edges, dims = 1)...)

ip2 = inpoly2(permutedims(traj2d.x0), nodes, edges) 
# stats[:, 1:2, area], where : has the index of the points and 1:2 are [inside, onboundary]
# therefore, stats[:, 1, k] is a bitvector of point membership to the kth polygon