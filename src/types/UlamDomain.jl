"""
    UlamDomain{S, T, U}

The domain defining the Ulam problem setup.

### Fields
- `corners`: The rectangular bounding box of the polygon in `domain`. Binning algorithms start by covering this rectangle.
- `domain`: Data outside this are considered in nirvana.
- `poly_type`: The kind of polygons which cover `domain`.
- `poly_number`: The number of polygons which cover `domain`.
- `stoc_type`: The manner in which trajectories which leave `domain` are reinjected.
- `stoc_source`: The location in which to reinject data for the `"source"` reinjection algorithm.
"""
struct UlamDomain{S<:AbstractString, T<:Real, U<:Integer}
    corners::Vector{T}
    domain::UlamPolygon{T}
    poly_type::S
    poly_number::U
    stoc_type::S
    stoc_source::Union{UlamPolygon{T}, Matrix{T}, Nothing}
    rseed::U
end


"""
    UlamDomain(domain; [...])

Construct an `UlamDomain` defined by the `UlamPolygon` in `domain`.

Points outside `domain` will be considered in nirvana.

### Optional Arguments
- `poly_type`: One of `"rec"`, `"sqr"`, `"tri"`, `"hex"`, and `"vor"` for coverings by rectangles, squares, hexagons or Voronoi tesselation. The default is rectangles.
- `poly_number`: The number of polygons requested. The default is 500 for rectangles/squares/hexagons and 100 for Voronoi.
- `stoc_type`: Picks the stochasticization algorithm; one of `"data"` or `"source"`. The default is data.
- `stoc_source`::`UlamPolygon`: Polygons in the covering that intersect `stoc_source` will have data re-injected uniformly through them in the `source` algorithm. 
- `stoc_source`::`Matrix`: Polygons in the covering which contain the points in the `stoc_source` matrix will have data re-injected uniformly through them in the `source` algorithm. 
- `rseed`: A seed for reproducing the random initialization of the kmeans algorithm in the Voronoi covering.
"""
function UlamDomain(
    domain::UlamPolygon{T};
    poly_type::S = global_poly_types[1],
    poly_number::U = global_poly_number_default[poly_type],
    stoc_type::S = global_stoc_types[1],
    stoc_source::Union{UlamPolygon{T}, Matrix{T}, Nothing} = nothing,
    rseed::U = global_rseed_default) where {S<:AbstractString, T<:Real, U<:Integer}

    xmin, xmax = extrema(domain.nodes[:, 1])
    ymin, ymax = extrema(domain.nodes[:, 2])

    @assert xmax > xmin
    @assert ymax > ymin

    corners::Vector{Float64} = [xmin, xmax, ymin, ymax]

    @assert poly_type in global_poly_types
    @assert poly_number > 1
    @assert stoc_type in global_stoc_types

    if typeof(stoc_source) <: AbstractMatrix
        @assert size(stoc_source, 1) > 0 
        @assert size(stoc_source, 2) == 2
    end

    if stoc_source !== nothing && stoc_type != "source"
        # assume the user wants to use the source algorithm if they provide a stoc_source even if they
        # didn't specify "source"
        stoc_type = "source"
    elseif stoc_source === nothing && stoc_type == "source"
        @warn "The `source`` reinjection algorithm was requested, but `stoc_source` was not provided. The 
        `stoc_source`` will be set to the domain. This is equivalent to reinjecting data uniformly across all boxes."
        stoc_source = domain
    end

    return UlamDomain{String, Float64, Int64}(
        corners,
        domain,
        poly_type,
        poly_number,
        stoc_type,
        stoc_source,
        rseed)
end

"""
    UlamDomain(xmin, xmax, ymin, ymax; [...])

Construct an `UlamDomain` defined by the rectangle with bottom left corner (`xmin`, `ymin`) and top right corner (`xmax`, `ymax`).

This is equivalent to `UlamDomain(domain; [...])` where the provided domain is a rectangular `UlamPolygon`.
"""
function UlamDomain(
    xmin::T,
    xmax::T,
    ymin::T,
    ymax::T;
    poly_type::S = global_poly_types[1],
    poly_number::U = global_poly_number_default[poly_type],
    stoc_type::S = global_stoc_types[1],
    stoc_source::Union{UlamPolygon{T}, Matrix{T}, Nothing} = nothing,
    rseed::U = global_rseed_default) where {S<:AbstractString, T<:Real, U<:Integer}

    @assert xmax > xmin
    @assert ymax > ymin

    corners::Vector{Float64} = [xmin, xmax, ymin, ymax]
    domain = UlamPolygon([xmin ymin; xmin ymax; xmax ymax; xmax ymin])

    return UlamDomain{String, Float64, Int64}(
        corners,
        domain,
        poly_type,
        poly_number,
        stoc_type,
        stoc_source,
        rseed)
end

function Base.show(io::IO, x::UlamDomain)
    print(io, "UlamDomain[")
    show(io, x.corners)
    print(io, ", ")
    show(io, x.poly_number)
    print(io, "@")
    show(io, x.poly_type)
    if x.domain !== nothing
        print(io, " w/ domain")
    end
    print("]")
end