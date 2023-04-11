using Graphs:SimpleDiGraph,  strongly_connected_components

"""
    ulam_nirvana(traj, domain, polys)

Execute the core algorithm of Ulam's method with a nirvana state. Return an [`UlamResult`](@ref).
"""
function ulam_nirvana(traj::UlamTrajectories, domain::UlamDomain, polys::Vector{UlamPolygon{T}}) where {T<:Real}
    ############ Assign polygon indices to trajectory points

    # find the polygon indices of initial and final trajectory points 
    ip0 = inpoly([traj.x0 ;; traj.y0], PolyTable(polys));
    ipT = inpoly([traj.xT ;; traj.yT], PolyTable(polys));

    # discard polygons that don't contain data
    polys = polys[[i for i in keys(ip0.contains)]]
    n_polys = length(polys)

    ############ Construct P_closed

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

    counts = [Int64(sum(P_closed[i, :])) for i = 1:n_polys]

    ############ Find the largest strongly connected component

    # create the adjacency matrix of P_closed; note that nirvana is excluded since we assume it's always connected
    Padj::Matrix{Int64} = [P_closed[i,j] != 0.0 ? 1 : 0 for i in 1:n_polys, j in 1:n_polys]

    # Construct the directed graph with adjacency matrix Padj and find its sccs; sort to get the largest scc.
    scc = sort(strongly_connected_components(SimpleDiGraph(Padj)), by = length)[end]

    # In general, scc itself is not sorted
    # we sort one more time and put nirvana back at the end to get the largest scc of P_closed
    # P_closed, polys and n_polys are updated to refer to the scc; store the disconnected polys in polys_dis
    largest_component = sort(scc)
    P_closed = P_closed[[largest_component; n_polys + 1], [largest_component; n_polys + 1]]
    n_polys = length(largest_component)
    polys_dis::Vector{eltype(polys)} = [polys[i] for i=1:n_polys if !(i in largest_component)]
    polys = polys[largest_component]

    ############ stochasticizing/reinjection

    # the "data" algorithm is applied by default and needs no further calculations

    # now we handle the reinjection if the user requested source
    if domain.stoc_type == "source" && sum(P_closed[n_polys + 1, :]) > 0.0
        # findind the polygons that intersect with the given sto_polygon
        in_source = [ulam_intersects(domain.sto_polygon, polys[i]) for i in largest_component]

        # all the reinjection counts are redistributed evenly at those locations
        total_counts = sum(P_closed[n_polys + 1, :])
        P_closed[n_polys + 1, :] = zeros(n_polys + 1)
        P_closed[n_polys + 1, :][in_source] .= total_counts/length(in_source)
    end

    # stochasticizing P_closed, except last row which is handled separately next
    for i = 1:n_polys
        P_closed[i, :] = P_closed[i, :]/sum(P_closed[i, :])
    end

    # reinjection trajectories
    if sum(P_closed[n_polys + 1, :]) > 0.0  
        P_closed[n_polys + 1, :] = P_closed[n_polys + 1, :]/sum(P_closed[n_polys + 1, :])
    else
        @info "There are no trajectories from nirvana to the interior."
    end

    if sum(P_closed[:, n_polys + 1]) == 0.0
        @info "There are no trajectories from the interior to nirvana."
    end

    ############ collect results and info, and return

    n_polys_requested = domain.poly_number
    n_polys_no_data = n_polys_requested - length(keys(ip0.contains))
    n_polys_dis = length(polys_dis)
    n_counts_total = sum(counts) 
    n_counts_removed_scc = n_counts_total - sum(counts[largest_component])
    poly_type = domain.poly_type
    info = UlamInfo(n_polys_requested, n_polys_no_data, n_polys_dis, n_counts_total, n_counts_removed_scc, poly_type)

    counts = counts[largest_component]

    return UlamResult(P_closed, polys, polys_dis, counts, info)
end