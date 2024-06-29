"""
    struct UlamMatrix
"""
struct UlamMatrix
    O2O::Matrix{Float64}
    O2ω::Vector{Float64}
    ω2O::Vector{Float64}
end

P_open(P::UlamMatrix) = P.O2O

P_closed(P::UlamMatrix) = vcat(hcat(P.O2O, P.O2ω), [P.ω2O ; 0.0]')