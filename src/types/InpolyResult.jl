struct InpolyResult{U<:Integer}
    inds::Vector{U}
    contains::Dict{U, U}
end