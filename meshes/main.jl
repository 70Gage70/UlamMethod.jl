using Meshes, CairoMakie
using ArgCheck
using PolygonInbounds
using Graphs:SimpleDiGraph,  strongly_connected_components

using Distributions # FOR TESTING ONLY

include("boundary.jl")
include("bins.jl")
include("traj.jl")
include("membership.jl")
include("reinjection.jl")

### 1d
boundary1d = Segment((0,), (1,)) |> Boundary
bin1d = bin(boundary1d, LineBinner(500))
traj1d = Trajectories(rand(1, 10), rand(1, 10))
membs1d = membership(traj1d, bin1d)

### 2d
boundary2d = Ngon([(0,0),(6,0),(1,7),(1,6)]...) |> Boundary
bin2d = bin(boundary2d, RectangleBinner(500))
x0_rand = rand(MvNormal([1, 4], [1 0; 0 1]), 1000)
xT_rand = x0_rand + rand(MvNormal([0, 0], [1 0; 0 1]), size(x0_rand, 2))
traj2d = Trajectories(x0_rand, xT_rand)
membs2d = membership(traj2d, bin2d)

n_bins = length(bin2d.bins)
n_points = size(traj2d.x0, 2)
Pij = zeros(n_bins + 1, n_bins + 1)
x0_idx, xT_idx = membs2d

### COMPUTE TRANSITION MATRIX
for i = 1:n_points
    x0, xT = membs2d[1][i], membs2d[2][i]
    if !isnothing(x0) && !isnothing(xT) # interior to interior
        Pij[x0, xT] += 1
    elseif !isnothing(x0) && isnothing(xT) # interior to nirvana
        Pij[x0, end] += 1
    elseif isnothing(x0) && !isnothing(xT) # nirvana to interior
        Pij[end, xT] += 1
    end # otherwise, transition from nirvana to nirvana and ignore
end 

### LARGEST STRONGLY CONNECTED COMPONENT 

# create the adjacency matrix of Pij; note that nirvana is excluded since we assume it's always connected
Padj = [iszero(Pij[i,j]) ? 0 : 1 for i in 1:n_bins, j in 1:n_bins]

scc = strongly_connected_components(SimpleDiGraph(Padj)) # Construct the directed graph with adjacency matrix Padj
scc = sort(scc, by = length)[end] # get the largest scc
scc = [sort(scc) ; n_bins + 1] # ensure bin labels are sorted, then put nirvana back

Pij = Pij[scc, scc]
bins_full = Bins(splice!(bin2d.bins, scc[1:end-1]))
bins_empty = bin2d

### REINJECTION ALGORITHM
Pij = reinject(bin2d, Pij, DataReinjection())

### STOCHASTICIZE
Pij = Pij ./ sum(Pij, dims = 2)
