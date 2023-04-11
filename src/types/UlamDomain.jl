"""
    UlamDomain{S, T, U}

The domain defining the Ulam problem setup.
"""
struct UlamDomain{S<:AbstractString, T<:Real, U<:Integer}
    corners::Vector{T}
    domain::Union{UlamPolygon{T}, Nothing}
    poly_type::S
    poly_number::U
    stoc_type::S
    stoc_polygon::Union{UlamPolygon{T},Nothing}
    rseed::U
end


"""
    UlamDomain(xmin, xmax, ymin, ymax; [...])

Construct an `UlamDomain` defined by the rectangle with bottom left corner (`xmin`, `ymin`) and top right corner (`xmax`, `ymax`).

### Optional Arguments
- `domain`: Points inside the rectangle, but outside the outside `domain` will be considered in nirvana. This can be used to refine the shape of the computational domain to an arbitrary `UlamPolygon`, not just the default rectangle.
- `poly_type`: One of `"sqr"`, `"hex"`, and `"vor"` for coverings by squares, hexagons or Voronoi tesselation. The default is squares.
- `poly_number`: The number of polygons requested. The default is 500 for squares/hexagons and 100 for Voronoi.
- `stoc_type`: Picks the stochasticization algorithm; one of `"data"` or `"source"`. The default is data.
- `stoc_polygon`: Polygons in the covering that intersect `stoc_polygon` will have data re-injected uniformly through them in the `source` algorithm. 
- `rseed`: A seed for reproducing the random initialization of the kmeans algorithm in the Voronoi covering.
"""
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
    rseed::U = global_rseed_default) where {S<:AbstractString, T<:Real, U<:Integer}

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

    if stoc_polygon !== nothing && stoc_type != "source"
        # assume the user wants to use the source algorithm if they provide a stoc_polygon even if they
        # didn't specify "source"
        stoc_type = "source"
    elseif stoc_polygon === nothing && stoc_type == "source"
        @warn "The `source`` reinjection algorithm was requested, but `stoc_polygon` was not provided. The 
        `stoc_polygon`` will be set to a rectangle with edges defined by `corners`. This is equivalent to
        reinjecting data uniformly across all boxes."
        stoc_polygon = UlamPolygon([xmin ymin; xmin ymax; xmax ymax; xmin ymax])
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
    if x.domain !== nothing
        print(io, " w/ domain")
    end
    print("]")
end