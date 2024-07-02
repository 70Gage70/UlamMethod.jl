"""
    ulam_method(traj, boundary, binner; reinj_algo)

Run the main Ulam's method calculation and return an [`UlamResult`](@ref).

### Arguments

- `traj`: A [`Trajectories`](@ref) object, holding the short-range trajectory data.
- `binner`: A [`BinningAlgorithm`](@ref) that specifies the algorithm used to partition the boundary into bins.

### Optional Arguments

- `reinj_algo`: A [`ReinjectionAlgorithm`](@ref) that specifies how trajectories pointing \
from nirvana to the interior should be reinjected. Default [`DataReinjection`](@ref).
"""
function ulam_method(
    traj::Trajectories{Dim},
    binner::BinningAlgorithm{Dim};
    reinj_algo::ReinjectionAlgorithm = DataReinjection()) where {Dim}

    ### COMPUTE BINS BASED ON BOUNDARY
    bins = binner.bins.bins
    n_bins = length(bins)
    n_bins_initial = n_bins
    pos2idx = collect(1:n_bins)
    
    ### COMPUTE BIN MEMBERSHIP
    x0_idx, xT_idx = membership(traj, binner)

    ### COMPUTE TRANSITION MATRIX
    Pij = zeros(n_bins + 1, n_bins + 1)
    n_points = size(traj.x0, 2)
    for i = 1:n_points
        x0, xT = x0_idx[i],xT_idx[i]
        if !isnothing(x0) && !isnothing(xT) # interior to interior
            Pij[x0, xT] += 1
        elseif !isnothing(x0) && isnothing(xT) # interior to nirvana
            Pij[x0, end] += 1
        elseif isnothing(x0) && !isnothing(xT) # nirvana to interior
            Pij[end, xT] += 1
        end # otherwise, transition from nirvana to nirvana and ignore
    end 

    ### REMOVE EMPTY BINS, POS2IDX
    full = [findall(!iszero, vec(sum(Pij[1:n_bins, 1:n_bins], dims = 2))) ; n_bins + 1] # don't check nirvana
    not_full = setdiff(1:n_bins+1, full)
    Pij = Pij[full, full]
    splice!(bins, not_full)
    n_bins = length(bins)

    pos2idx = pos2idx[full[1:end-1]] # post2idx[i] gives the original index of what is now the ith bin
    
    ### LARGEST STRONGLY CONNECTED COMPONENT, POS2IDX

    # create the adjacency matrix of Pij; note that nirvana is excluded since we assume it's always connected
    Padj = [iszero(Pij[i,j]) ? 0 : 1 for i in 1:n_bins, j in 1:n_bins]

    scc = strongly_connected_components(SimpleDiGraph(Padj)) # Construct the directed graph with adjacency matrix Padj
    scc = sort(scc, by = length)[end] # get the largest scc
    scc = [sort(scc) ; n_bins + 1] # ensure bin labels are sorted, then put nirvana back
    not_scc = setdiff(1:n_bins+1, scc)

    Pij = Pij[scc, scc]
    bins_dis = splice!(bins, not_scc)
    n_bins = length(bins)

    pos2idx = pos2idx[scc[1:end-1]]
    binner.idx2pos .= [findfirst(x -> x == i, pos2idx) for i = 1:n_bins_initial]

    ### VALIDATION
    if any(iszero.(vec(sum(Pij[1:n_bins, 1:n_bins], dims = 2))))
        error("The transition probability matrix contains rows with no counts. This probably means that the trajectories do not create enough communication between bins.")
    end

    ### REINJECTION ALGORITHM
    reinject!(binner, Pij, reinj_algo)

    ### STOCHASTICIZE
    Pω2O = Pij[n_bins+1, 1:n_bins]
    try
        Pω2O = Pω2O / sum(Pω2O)
    catch
        # Pω2O is a vector of zeros, so leave it alone
    end

    Pij = Pij[1:end-1, :] ./ sum(Pij[1:end-1, :], dims = 2)
    PO2O = Pij[1:n_bins, 1:n_bins]
    PO2ω = Pij[1:n_bins, n_bins+1]

    ### RETURN
    return UlamResult(PO2O, PO2ω, Pω2O, binner, Bins(bins_dis))
end