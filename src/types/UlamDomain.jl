struct UlamDomain{S<:AbstractString, T<:Real, U<:Integer}
    corners::Vector{T}
    domain::Union{UlamPolygon{T}, Nothing}
    poly_type::S
    poly_number::U
    stoc_type::S
    stoc_polygon::Union{UlamPolygon{T},Nothing}
    rseed::U
end

function UlamDomain(
    xmin::T,
    xmax::T,
    ymin::T,
    ymax::T;
    domain::Union{UlamPolygon{T}, Nothing} = nothing,
    poly_type::S = global_poly_types[1],
    poly_number::U = global_poly_number_default[poly_type],
    stoc_type::S = global_stoc_types[1],
    stoc_polygon::Union{UlamPolygon{T}, Nothing} = nothing,
    rseed::U = 123) where {S<:AbstractString, T<:Real, U<:Integer}

    @assert xmax > xmin
    @assert ymax > ymin

    corners::Vector{Float64} = [xmin, xmax, ymin, ymax]

    # if corners is ever replaced by a generic domain
    # if domain === nothing
    #     domain = UlamPolygon([xmin ymin; xmin ymax; xmax ymax; xmin ymax])
    # end

    @assert poly_type in global_poly_types
    @assert poly_number > 1
    @assert stoc_type in global_stoc_types

    if stoc_polygon != nothing && stoc_type != "source"
        # assume the user wants to use the source algorithm if they provide a stoc_polygon even if they
        # didn't specify "source"
        stoc_type = "source"
    end

    return UlamDomain{String, Float64, Int64}(
        corners,
        domain,
        poly_type,
        poly_number,
        stoc_type,
        stoc_polygon,
        rseed)
end

function Base.show(io::IO, x::UlamDomain)
    print(io, "UlamDomain[")
    show(io, x.corners)
    print(io, ", ")
    show(io, x.poly_number)
    print(io, "@")
    show(io, x.poly_type)
    if x.domain != nothing
        print(io, " w/ domain")
    end
    print("]")
end