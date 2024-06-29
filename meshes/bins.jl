"""
    struct Bins{K, Dim, CRS}

A container type for bins.

### Fields

- `bins`: A vector of `Polytope` objects of dimension `Dim` and coordinate reference system `CRS`.
"""
struct Bins{K, Dim, CRS}
    bins::Vector{<:Polytope{K, Dim, CRS}}
end

"""
    abstract type BinningAlgorithm{Dim}

An abstract type for binning algorithms of dimension `Dim`.

Each subtype `binner` should implement the following:

- a function `_bin(boundary, binner)` that partitions `boundary` and return `Bins`.
"""
abstract type BinningAlgorithm{Dim} end

### 1D ALGORITHMS

"""
    struct LineBinner

Bin a one dimensional `Segment` (line segement) with `nbins` equally-spaced bins.

### Fields

- `nbins`: The number of bins.
- `hardclip`: When set, bins are clipped to the boundary so there is no overlap.

### Constructor 

`LineBinner(nbins; hardclip = true)`
"""
struct LineBinner <: BinningAlgorithm{1}
    nbins::Int64
    hardclip::Bool

    function LineBinner(nbins::Int64; hardclip::Bool = true)
        return new(nbins, hardclip)
    end
end

function _bin(boundary::Boundary{K, 1, CRS}, binner::LineBinner) where {K, CRS}
    nbins = binner.nbins
    bbox = Meshes.boundingbox(boundary.boundary)
    grid = CartesianGrid(bbox.min, bbox.max, dims = (nbins, ))

    bins = Polytope{K, 1, CRS}[]
    for bin_ in elements(grid)
        isect = binner.hardclip ? intersect(bin_, boundary.boundary) : intersects(bin_, boundary.boundary) ? bin_ : nothing
        if isect isa Segment
            push!(bins, isect)
        end
    end

    @info bins

    return Bins(bins)
end 

### 2D ALGORITHMS

"""
    struct RectangleBinner

Bin a two dimensional `Polygon` with a tight covering of rectangles. The rectangles
are chosen to be as close as possible to squares and so the final number of bins
may be slightly different than the number requested.

### Fields

- `nbins`: The number of bins requested.
- `hardclip`: When set, bins are clipped to the boundary so there is no overlap.

### Constructor 

`RectangleBinner(nbins; hardclip = true)`
"""
struct RectangleBinner <: BinningAlgorithm{2}
    nbins::Int64
    hardclip::Bool

    function RectangleBinner(nbins::Int64; hardclip::Bool = true)
        return new(nbins, hardclip)
    end
end

function _bin(boundary::Boundary{K, 2, CRS}, binner::RectangleBinner) where {K, CRS}
    nbins = binner.nbins

    bbox = Meshes.boundingbox(boundary.boundary)
    bbox_W = coords(bbox.max).x - coords(bbox.min).x
    bbox_L = coords(bbox.max).y - coords(bbox.min).y
    
    n_x = sqrt(bbox_W*nbins/bbox_L) |> x -> round(Int64, x)
    n_y = sqrt(bbox_L*nbins/bbox_W) |> x -> round(Int64, x)
    
    grid = CartesianGrid(bbox.min, bbox.max, dims = (n_x, n_y))
    bins = Polytope{K, 2, CRS}[]

    for bin_ in elements(grid)
        isect = binner.hardclip ? intersect(bin_, boundary.boundary) : intersects(bin_, boundary.boundary) ? bin_ : nothing
        if isect isa PolyArea
            # if the intersection is a PolyArea (i.e. a polygon possibly with holes), we convert it to an Ngon with the outer ring
            push!(bins, Ngon(isect.rings[1].vertices...))
        elseif isect isa Ngon
            push!(bins, isect)
        end
    end

    return Bins(bins)
end 

"""
    bin(boundary, binner)

Bin the `boundary` according to the `binner` algorithm.
"""
function bin(boundary::Boundary{K, Dim, CRS}, binner::BinningAlgorithm{Dim}) where {K, Dim, CRS}
    return _bin(boundary, binner)
end