struct UlamResult{S<:AbstractString, T<:Real, U<:Integer}
    P_closed::Matrix{T}
    P_open::SubArray{T, 2, Matrix{T}, Tuple{UnitRange{U}, UnitRange{U}}, false}
    pi_closed::Vector{T}
    pi_open::SubArray{T, 1, Vector{T}, Tuple{UnitRange{U}}, true}
    polys::Vector{UlamPolygon{T}}
    polys_dis::Vector{UlamPolygon{T}}
    counts::Vector{U}
    info::UlamInfo{S, U}
end


function UlamResult(
    P_closed::Matrix{T},
    polys::Vector{UlamPolygon{T}},
    polys_dis::Vector{UlamPolygon{T}},
    counts::Vector{U},
    info::UlamInfo{S, U}) where {S<:AbstractString, T<:Real, U<:Integer}

    @assert size(P_closed, 1) == size(P_closed, 2)
    @assert size(P_closed, 1) > 0

    @assert length(polys) > 0
    
    pi_closed = abs.(LinearAlgebra.normalize(LinearAlgebra.eigvecs(transpose(P_closed))[:,size(P_closed)[1]], 1))

    P_open = view(P_closed, 1:size(P_closed, 1) - 1, 1:size(P_closed, 1) - 1)
    pi_open = view(pi_closed, 1:size(pi_closed, 1) - 1)

    return UlamResult{S, T, U}(P_closed, P_open, pi_closed, pi_open, polys, polys_dis, counts, info)
end


function Base.show(io::IO, x::UlamResult)
    print(io, " UlamResult")
    println(io)
    print("  Polys requested: ")
    show(io, x.info.n_polys_requested)
    println(io)
    print("  Polys type: ")
    show(io, x.info.poly_type)
    println(io)
    print("  Polys obtained: ")
    show(io, length(x.polys))
    println(io)
    print("  Polys disconnected: ")
    show(io, length(x.polys_dis))
end