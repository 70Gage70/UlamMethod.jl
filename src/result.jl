"""
    struct UlamResult{K, Dim, CRS}

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
`bins`: The bins corresponding to the transition probability matrix.
`bins_dis`: The bins that contained data but were removed when the largest strongly connected component was taken.

### Methods

Use `P_open(UlamResult)` to access `P_O2O` and `P_closed(UlamResult)` to access the full matrix.

Use `bins(UlamResult)` and `bins_dis(UlamResult)` to access `bins` and `bins_dis`, respectively.
"""
struct UlamResult{K, Dim, CRS}
    P_O2O::Matrix{Float64}
    P_O2ω::Vector{Float64}
    P_ω2O::Vector{Float64}
    bins::Bins{K, Dim, CRS}
    bins_dis::Bins{K, Dim, CRS}
end

"""
    P_open(ur)

Return `ur.O2O`, i.e. the transition matrix without nirvana.
"""
P_open(ur::UlamResult) = ur.P_O2O

"""
    P_closed(ur)

Return the full transition matrix.
"""
P_closed(ur::UlamResult) = vcat(hcat(ur.P_O2O, ur.P_O2ω), [ur.P_ω2O ; 0.0]')

"""
    bins(ur)

Return the bins associated to the transition matrix.
"""
bins(ur::UlamResult) = ur.bins

"""
    bins_dis(ur)

Return the bins that contained data but were removed when the largest strongly connected component was taken.
"""
bins_dis(ur::UlamResult) = ur.bins_dis