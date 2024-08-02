"""
    struct UlamResult{Dim, M, CRS}

A container for the result of the main Ulam's method calculation. 

This consists of three main components, the transition probability matrix, the associated 
bins and the bins that were disconnected.

The transition probability matrix is of the form

| `P_O2O` | `P_O2ω` | 

| `P_ω2O` |  `0`  |

Where `O` represents the "open" (interior) region and `ω` reprsents the "nirvana" (exterior) region.
The union of `O` and `ω` form a closed region.

### Fields 

- `P_O2O`: As in diagram.
- `P_O2ω`: As in diagram.
- `P_ω2O`: As in diagram.
- `binner`: The [`BinningAlgorithm`](@ref) used to bin the data.
- `bins_dis`: The bins that contained data but were removed when the largest strongly connected component was taken.

### Methods

    P_closed(UlamResult)

Access the full matrix.

    P_open(UlamResult)
    
Access `P_O2O`.

    bins(UlamResult)
    
Access `bins`.

    bins_dis(UlamResult)

Access `bins_dis`.

    membership(data, ulam_result)

Compute `membership(data, ulam_result.binner)`.
"""
struct UlamResult{Dim, M, CRS}
    P_O2O::Matrix{Float64}
    P_O2ω::Vector{Float64}
    P_ω2O::Vector{Float64}
    binner::BinningAlgorithm{Dim}
    bins_dis::Bins{Dim, M, CRS}
end

"""
    P_open(ulam_result)

Return `ulam_result.O2O`, i.e. the transition matrix without nirvana.
"""
P_open(ur::UlamResult) = ur.P_O2O

"""
    P_closed(ulam_result)

Return the full transition matrix.
"""
P_closed(ur::UlamResult) = vcat(hcat(ur.P_O2O, ur.P_O2ω), [ur.P_ω2O ; 0.0]')

"""
    bins(ulam_result)

Return the bins associated to the transition matrix.
"""
bins(ur::UlamResult) = ur.binner.bins

"""
    bins_dis(ulam_result)

Return the bins that contained data but were removed when the largest strongly connected component was taken.
"""
bins_dis(ur::UlamResult) = ur.bins_dis

"""
    membership(data, binner)

Takes a `Dim x N_points` matrix of points and returns a vector `memb` where `memb[i] = j` if `data[:,i]` is 
inside `binner.bins[j]` and `memb[i] = nothing` if `data[:,i]` is not inside any bin.

    membership(traj, binner)

Compute `(membership(traj.x0, binner), membership(traj.xT, binner))`.

    membership(data, ulam_result)

Compute `membership(data, ulam_result.binner)`.
"""
function membership(data::Matrix{<:Real}, ur::UlamResult{Dim, M, CRS}) where {Dim, M, CRS}
    @argcheck size(data, 1) == Dim
    return membership(data, ur.binner)
end

