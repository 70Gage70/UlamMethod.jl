"""
    struct Bins{Dim, CRS}

A container type for bins.

### Fields

- `bins`: A vector of `Geometry` objects of dimension `Dim` and coordinate reference system `CRS`.
"""
struct Bins{Dim, CRS}
    bins::Vector{<:Geometry{Dim, CRS}}
end

"""
    abstract type AbstractBinner{Dim}

An abstract type for binning algorithms of dimension `Dim`.

Each subtype `binner` should implement a function `_bin(boundary, binner)` that
returns `Bins`.
"""
abstract type AbstractBinner{Dim} end

### 1D ALGORITHMS

"""
    struct LineBinner

Bin a one dimensional `Chain` (line segement) with `nbins` equally-spaced bins.

### Fields

- `nbins`: The number of bins.
"""
struct LineBinner <: AbstractBinner{1}
    nbins::Int64
end

function _bin(boundary::Chain, binner::LineBinner)
    nbins = binner.nbins
    bbox = Meshes.boundingbox(boundary)
    grid = CartesianGrid(bbox.min, bbox.max, dims = (nbins, ))
    return Bins([elem for elem in elements(grid) if intersects(elem, boundary)])
end 

### 2D ALGORITHMS

"""
    struct RectangleBinner

Bin a two dimensional `Polygon` with a tight covering of rectangles. The rectangles
are chosen to be as close as possible to squares and so the final number of bins
may be slightly different than the number requested.

### Fields

- `nbins`: The number of bins requested.
"""
struct RectangleBinner <: AbstractBinner{2}
    nbins::Int64
end

function _bin(boundary::Polygon, binner::RectangleBinner)
    nbins = binner.nbins

    bbox = Meshes.boundingbox(boundary)
    bbox_W = coords(bbox.max).x - coords(bbox.min).x
    bbox_L = coords(bbox.max).y - coords(bbox.min).y
    
    n_x = sqrt(bbox_W*nbins/bbox_L) |> x -> round(Int64, x)
    n_y = sqrt(bbox_L*nbins/bbox_W) |> x -> round(Int64, x)
    
    grid = CartesianGrid(bbox.min, bbox.max, dims = (n_x, n_y))
    return Bins([elem for elem in elements(grid) if intersects(elem, boundary)])
end 

"""
    bin(boundary, binner; hardclip)

Bin the `boundary` according to the `binner` algorithm.

### Optional Arguments

- `hardclip`: When set, clip all bins to the boundary so there is no overlap. Default: `false`.
"""
function bin(
    boundary::Polytope{K, Dim, CRS}, 
    binner::AbstractBinner{K};
    hardclip::Bool = false) where {K, Dim, CRS}

    bins = _bin(boundary, binner)
    return hardclip ? Bins([intersect(boundary, bin) for bin in bins]) : bins
end