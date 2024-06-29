"""
    ulam_method(traj, boundary, binner; reinj_algo)

Run the main Ulam's method calculation and return an [`UlamResult`](@ref).

### Arguments

- `traj`: A [`Trajectories`](@ref) object, holding the short-range trajectory data.
- `boundary`: A [`Boundary`](@ref) object, holding the geometry that defines the computational boundary.
- `binner`: A [`BinningAlgorithm`](@ref) that specifies the algorithm used to partition the boundary into bins.

### Optional Arguments

- `reinj_algo`: A [`ReinjectionAlgorithm`](@ref) that specifies how trajectories pointing \
from nirvana to the interior should be reinjected. Default [`DataReinjection`](@ref).
"""
function ulam_method(
    traj::Trajectories{Dim},
    boundary::Boundary{K, Dim, CRS}, 
    binner::BinningAlgorithm{Dim};
    reinj_algo::ReinjectionAlgorithm = DataReinjection()) where {K, Dim, CRS}

    ### COMPUTE BINS BASED ON BOUNDARY
    bins = bin(boundary, binner)
    n_bins = length(bins.bins)
    
    ### COMPUTE BIN MEMBERSHIP
    x0_idx, xT_idx = membership(traj, bins)

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

    ### REMOVE EMPTY BINS
    full = [findall(!iszero, vec(sum(Pij[1:n_bins, 1:n_bins], dims = 2))) ; n_bins + 1] # don't check nirvana
    Pij = Pij[full, full]
    bins_full = Bins(splice!(bins.bins, full[1:end-1]))
    n_bins = length(bins_full.bins)

    ### LARGEST STRONGLY CONNECTED COMPONENT 

    # create the adjacency matrix of Pij; note that nirvana is excluded since we assume it's always connected
    Padj = [iszero(Pij[i,j]) ? 0 : 1 for i in 1:n_bins, j in 1:n_bins]

    scc = strongly_connected_components(SimpleDiGraph(Padj)) # Construct the directed graph with adjacency matrix Padj
    scc = sort(scc, by = length)[end] # get the largest scc
    scc = [sort(scc) ; n_bins + 1] # ensure bin labels are sorted, then put nirvana back

    Pij = Pij[scc, scc]
    bins_final = Bins(splice!(bins_full.bins, scc[1:end-1]))
    bins_dis = bins_full
    n_bins = length(bins_final.bins)

    ### REINJECTION ALGORITHM
    Pij = reinject(bins_final, Pij, reinj_algo)

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
    return UlamResult(PO2O, PO2ω, Pω2O, bins_final, bins_dis)
end