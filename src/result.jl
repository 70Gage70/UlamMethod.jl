"""
    struct UlamResult{Dim, CRS}

A container for the result of the main Ulam's method calculation. 

This consists of three main components, the transition probability matrix, the associated 
bins and the bins that were disconnected.

The transition probability matrix is of the form

| P_O2O | P_O2ω | 

| P_ω2O |  0  |

### Fields 

`P_O2O`: As in diagram.
`P_O2ω`: As in diagram.
`P_ω2O`: As in diagram.
`binner`: The [`BinningAlgorithm`] used to bin the data.
`bins_dis`: The bins that contained data but were removed when the largest strongly connected component was taken.

### Methods

Use `P_open(UlamResult)` to access `P_O2O` and `P_closed(UlamResult)` to access the full matrix.

Use `bins(UlamResult)` and `bins_dis(UlamResult)` to access `bins` and `bins_dis`, respectively.
"""
struct UlamResult{Dim, CRS}
    P_O2O::Matrix{Float64}
    P_O2ω::Vector{Float64}
    P_ω2O::Vector{Float64}
    binner::BinningAlgorithm{Dim}
    bins_dis::Bins{Dim, CRS}
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
    membership(data, ulam_result)

Takes a `Dim x N_points` matrix of points and returns a vector `memb` where `memb[i] = j` if `data[:,i]` is 
inside `ulam_result.binner.bins[j]` and `memb[i] = nothing` if `data[:,i]` is not inside any bin.
"""
function membership(data::Matrix{<:Real}, ur::UlamResult{Dim, CRS}) where {Dim, CRS}
    @argcheck size(data, 1) == Dim
    return membership(data, ur.binner)
end