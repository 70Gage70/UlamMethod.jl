"""
    InpolyResult{U}

Hold the result of [`inpoly``](@ref).
"""
struct InpolyResult{U<:Integer}
    inds::Vector{U}
    contains::Dict{U, U}
end