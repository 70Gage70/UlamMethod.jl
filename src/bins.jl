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
    bbox = Meshes.boundingbox(boundary.boundary)
    bbox_W = coords(bbox.max).x - coords(bbox.min).x
    bbox_L = coords(bbox.max).y - coords(bbox.min).y

    nbins = binner.hardclip ? ceil(Int64, binner.nbins*area(bbox)/area(boundary.boundary)) : binner.nbins
    
    n_x = sqrt(bbox_W*nbins/bbox_L) |> x -> round(Int64, x)
    n_y = sqrt(bbox_L*nbins/bbox_W) |> x -> round(Int64, x)
    
    grid = CartesianGrid(bbox.min, bbox.max, dims = (n_x, n_y))
    bins = Polytope{K, 2, CRS}[]

    for bin_ in elements(grid)
        isect = binner.hardclip ? intersect(bin_, boundary.boundary) : intersects(bin_, boundary.boundary) ? bin_ : nothing
        if isect isa PolyArea
            # if the intersection is a PolyArea (i.e. a polygon possibly with holes), we convert it to an Ngon with the outer ring
            verts = isect.rings[1].vertices
            if length(verts) >= 3
                push!(bins, Ngon(verts...))
            end
        elseif isect isa Ngon
            push!(bins, isect)
        end
    end

    return Bins(bins)
end 

"""
    struct HexagonBinner

Bin a two dimensional `Polygon` with a covering of regular hexagons. 

The final number of bins may be slightly different than the number requested.

### Fields

- `nbins`: The number of bins requested.
- `hardclip`: When set, bins are clipped to the boundary so there is no overlap.

### Constructor 

`HexagonBinner(nbins; hardclip = true)`
"""
struct HexagonBinner <: BinningAlgorithm{2}
    nbins::Int64
    hardclip::Bool

    function HexagonBinner(nbins::Int64; hardclip::Bool = true)
        return new(nbins, hardclip)
    end
end

function _bin(boundary::Boundary{K, 2, CRS}, binner::HexagonBinner) where {K, CRS}
    bbox = Meshes.boundingbox(boundary.boundary)
    xmin, xmax, ymin, ymax = coords(bbox.min).x.val, coords(bbox.max).x.val, coords(bbox.min).y.val, coords(bbox.max).y.val
    W = xmax - xmin
    L = ymax - ymin

    nbins = binner.hardclip ? ceil(Int64, binner.nbins*area(bbox)/area(boundary.boundary)) : binner.nbins
    
    α = L/W
    β = sqrt(1 + 24*sqrt(3)*α*nbins)
    n_x, n_y = ceil(Int64, (β - 1)/(4*sqrt(3)*α)), ceil(Int64, (β + 1)/6)
    R = max(W/(sqrt(3)*n_x), 2*L/(3*n_y - 1))
    r = sqrt(3)*R/2
    
    cs = []
    
    for y = 1:n_y
        x_left = isodd(y) ? 0.0 : -r
        n_x_corrected = isodd(y) ? n_x : n_x + 1
        for x = 1:n_x_corrected
            push!(cs, [x_left + (x - 1)*2*r, (y - 1)*3*R/2])
        end
    end 
    
    hexagons = []
    
    for c in cs
        push!(hexagons, [(c[1] + R*sin(π*k/3), c[2] + R*cos(π*k/3)) for k = 0:5])
    end
    
    hexagons = [Ngon(hex...) for hex in hexagons]
    
    c_box = ((xmin + xmax)/2, (ymin + ymax)/2)
    c_hex = centroid(boundingbox(hexagons)) |> p -> (coords(p).x.val, coords(p).y.val)
    hexagons = Translate(c_box .- c_hex)(hexagons)

    bins = Polytope{K, 2, CRS}[]
    for bin_ in hexagons
        isect = binner.hardclip ? intersect(bin_, boundary.boundary) : intersects(bin_, boundary.boundary) ? bin_ : nothing
        if isect isa PolyArea
            # if the intersection is a PolyArea (i.e. a polygon possibly with holes), we convert it to an Ngon with the outer ring
            verts = isect.rings[1].vertices
            if length(verts) >= 3
                push!(bins, Ngon(verts...))
            end
        elseif isect isa Ngon
            push!(bins, isect)
        end
    end

    return Bins(bins)
end 

"""
    struct TriangleBinner

Bin a two dimensional `Polygon` with a covering of equilateral triangles. 

The final number of bins may be slightly different than the number requested.

### Fields

- `nbins`: The number of bins requested.
- `hardclip`: When set, bins are clipped to the boundary so there is no overlap.

### Constructor 

`TriangleBinner(nbins; hardclip = true)`
"""
struct TriangleBinner <: BinningAlgorithm{2}
    nbins::Int64
    hardclip::Bool

    function TriangleBinner(nbins::Int64; hardclip::Bool = true)
        return new(nbins, hardclip)
    end
end

function _bin(boundary::Boundary{K, 2, CRS}, binner::TriangleBinner) where {K, CRS}
    bbox = Meshes.boundingbox(boundary.boundary)
    xmin, xmax, ymin, ymax = coords(bbox.min).x.val, coords(bbox.max).x.val, coords(bbox.min).y.val, coords(bbox.max).y.val
    W = xmax - xmin
    L = ymax - ymin

    nbins = binner.hardclip ? ceil(Int64, (1/6)*binner.nbins*area(bbox)/area(boundary.boundary)) : binner.nbins
    
    α = L/W
    β = sqrt(1 + 24*sqrt(3)*α*nbins)
    n_x, n_y = ceil(Int64, (β - 1)/(4*sqrt(3)*α)), ceil(Int64, (β + 1)/6)
    R = max(W/(sqrt(3)*n_x), 2*L/(3*n_y - 1))
    r = sqrt(3)*R/2
    
    cs = []
    
    for y = 1:n_y
        x_left = isodd(y) ? 0.0 : -r
        n_x_corrected = isodd(y) ? n_x : n_x + 1
        for x = 1:n_x_corrected
            push!(cs, [x_left + (x - 1)*2*r, (y - 1)*3*R/2])
        end
    end 
    
    triangles = []
    
    for c in cs, t = 1:6
        push!(triangles, [
            (c[1], c[2]), 
            (c[1] + R*sin(π*t/3), c[2] + R*cos(π*t/3)), 
            (c[1] + R*sin(π*(t + 1)/3), c[2] + R*cos(π*(t + 1)/3))
            ]
        )
    end
    
    triangles = [Ngon(hex...) for hex in triangles]
    
    c_box = ((xmin + xmax)/2, (ymin + ymax)/2)
    c_hex = centroid(boundingbox(triangles)) |> p -> (coords(p).x.val, coords(p).y.val)
    triangles = Translate(c_box .- c_hex)(triangles)

    bins = Polytope{K, 2, CRS}[]
    for bin_ in triangles
        isect = binner.hardclip ? intersect(bin_, boundary.boundary) : intersects(bin_, boundary.boundary) ? bin_ : nothing
        if isect isa PolyArea
            # if the intersection is a PolyArea (i.e. a polygon possibly with holes), we convert it to an Ngon with the outer ring
            verts = isect.rings[1].vertices
            if length(verts) >= 3
                push!(bins, Ngon(verts...))
            end
        elseif isect isa Ngon
            push!(bins, isect)
        end
    end

    return Bins(bins)
end 

### ND ALGORITHMS

"""
    struct VoronoiBinner{Dim}

Bin `Dim`-dimensional dataset based on a Voronoi tesellation of points generated by a k-means clustering of \
initial trajectory data (i.e. `Trajectories.x0`).

### Fields 

- `nbins`: The number of bins requested (number of k-means clusters).
- `traj`: The `Trajectories` data.
- `hardclip`: When set, bins are clipped to the boundary so there is no overlap.

### Constructor

`VoronoiBinner(nbins, traj; hardclip = true)`
"""
struct VoronoiBinner{Dim} <: BinningAlgorithm{Dim}
    nbins::Int64
    traj::Trajectories{Dim}
    hardclip::Bool
    
    function VoronoiBinner(nbins::Int64, traj::Trajectories{Dim}; hardclip::Bool = true) where {Dim}
        @argcheck Dim in [1, 2]
        new{Dim}(nbins, traj, hardclip)
    end
end

function _bin(boundary::Boundary{K, Dim, CRS}, binner::VoronoiBinner{Dim}) where {K, Dim, CRS}
    nbins = binner.nbins
    pts = [points(boundary) ;; binner.traj.x0]

    centers = kmeans(Yinyang(), pts, nbins).centers |> x -> [Point(x[:, i]...) for i = 1:size(x, 2)]

    if Dim == 1
        grid = RectilinearGrid([coords(c).x for c in centers])
    elseif Dim == 2
        grid = tesselate(centers, VoronoiTesselation())
    end

    bins = Polytope{K, 2, CRS}[]

    for bin_ in elements(grid)
        isect = binner.hardclip ? intersect(bin_, boundary.boundary) : intersects(bin_, boundary.boundary) ? bin_ : nothing
        if isect isa PolyArea
            # if the intersection is a PolyArea (i.e. a polygon possibly with holes), we convert it to an Ngon with the outer ring
            verts = isect.rings[1].vertices
            if length(verts) >= 3
                push!(bins, Ngon(verts...))
            end
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