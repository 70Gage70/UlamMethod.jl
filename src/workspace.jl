using .UlamTypes
import Graphs

include("ulam-binners.jl")
include("helpers.jl")

traj = UlamTrajectories("x0x5-NA-undrogued.mat");
domain = UlamDomain(-100, 15, -9, 39, poly_number = 500);
polys = square_binner(traj, domain);


function get_scc_inds(P_closed::Matrix{Float64})::Vector{Int64}
    Psize = length(P_closed[1,:])
    P_open_size = Psize - 1 # exclude nirvana when calculating scc, assumes one nirvana state in last entry
    
    Padj = Matrix{Int64}(undef, P_open_size, P_open_size) # the adjacency matrix defined by P
    for i = 1:P_open_size
        for j = 1:P_open_size
            if P_closed[i, j] == 0.0
                Padj[i,j] = 0 
            else
                Padj[i,j] = 1
            end
        end
    end
    
    Pgraph = SimpleDiGraph(Padj)
    scc = sort(strongly_connected_components(Pgraph), by = length) 
    largest_component = sort(scc[end])
    push!(largest_component, Psize) # put nirvana back at the end, should stay sorted

    return largest_component
end


# find the polygon indices of initial trajectory points 
ip0 = inpoly([traj.x0 ;; traj.y0], PolyTable(polys));

# discard polygons that don't contain data
contains_data = [i for i in keys(ip0.contains)]
polys = polys[contains_data]
n_polys = length(polys)

# find the polygon indices of final trajectory points
ipT = inpoly([traj.xT ;; traj.yT], PolyTable(polys));

# construct P_closed
P_closed = zeros(n_polys + 1, n_polys + 1) # one extra for nirvana

for i = 1:length(ip0.inds)
    i0, iT = ip0.inds[i], ipT.inds[i]
    if i0 != 0 # point starts in domain
        i0 = ip0.contains[i0] # i0 is still with respect to the uncleaned polygons, so we just fix that with the contains Dict
        if iT != 0 # point ends in domain
            P_closed[i0, iT] = P_closed[i0, iT] + 1
        elseif iT == 0 # point ends outside domain
            P_closed[i0, n_polys + 1] = P_closed[i0, n_polys + 1] + 1
        end
    elseif i0 == 0 # point starts outside domain
        if iT != 0 # point ends in domain (reinjection)
            P_closed[n_polys + 1, iT] = P_closed[n_polys + 1, iT] + 1
        end

        # else: point from nirvana to self; ignored 
    end
end

