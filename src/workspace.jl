using .UlamTypes
using Graphs:SimpleDiGraph,  strongly_connected_components

include("ulam-binners.jl")
include("helpers.jl")

traj = UlamTrajectories("x0x5-NA-undrogued.mat");
domain = UlamDomain(-100, 15, -9, 39, poly_number = 500);
polys = square_binner(traj, domain);

#########################################################################################################
#########################################################################################################
#########################################################################################################

### Assign polygon indices to trajectory points

# find the polygon indices of initial and final trajectory points 
ip0 = inpoly([traj.x0 ;; traj.y0], PolyTable(polys));
ipT = inpoly([traj.xT ;; traj.yT], PolyTable(polys));

# discard polygons that don't contain data
polys = polys[[i for i in keys(ip0.contains)]]
n_polys = length(polys)

### Construct P_closed

P_closed = zeros(n_polys + 1, n_polys + 1) # one extra for nirvana

for i = 1:length(ip0.inds)
    i0 = ip0.inds[i]
    iT = ipT.inds[i]

    if iT == 0 # (xT, yT) in nirvana
        if i0 != 0 # (x0, y0) in domain
            i0 = ip0.contains[i0] # index of cleaned polygon
            P_closed[i0, n_polys + 1] = P_closed[i0, n_polys + 1] + 1
        # elseif i0 == 0, (x0, y0) and (xT, yT) are both in nirvana; ignore
        end
    elseif iT in keys(ip0.contains) # (xT, yT) in the domain AND inside a box that contains data
        iT = ip0.contains[iT] 
        if i0 != 0 # (x0, y0) and (xT, yT) are both in domain
            i0 = ip0.contains[i0] 
            P_closed[i0, iT] = P_closed[i0, iT] + 1
        elseif i0 == 0 # (x0, y0) in nirvana but (xT, yT) in domain (reinjection)
            P_closed[n_polys + 1, iT] = P_closed[n_polys + 1, iT] + 1
        end        
    end
end

### Find the largest strongly connected component

# create the adjacency matrix of P_closed; note that nirvana is excluded since we assume it's always connected
Padj::Matrix{Int64} = [P_closed[i,j] != 0.0 ? 1 : 0 for i in 1:n_polys, j in 1:n_polys]

# Construct the directed graph with adjacency matrix Padj and find its sccs; sort to get the largest scc.
scc = sort(strongly_connected_components(SimpleDiGraph(Padj)), by = length)[end]

# In general, scc itself is not sorted, so we sort one more time and put nirvana back at the end to get the largest scc of P_closed
largest_component = [sort(scc); n_polys + 1]




P_closed = P_closed[largest_component, largest_component]

P_sto = [P_closed[i, j]/sum(P_closed[i,:]) for i in 1:n_polys, j in 1:n_polys];