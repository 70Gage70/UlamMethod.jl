"""
    struct Boundary{Dim, M, CRS}

A container type for the computational domain within which Ulam's method is applied.

The boundary has dimension `Dim` embedded in manifold `M` with coordinate reference system `CRS`.

The boundary is partitioned according to an [`BinningAlgorithm`](@ref).

### Fields

- `boundary`: A `Polytope` defining the boundary.

### Methods

    points(boundary)

Calculate a `Dim x N` matrix of boundary vertices.

### 1D Constructors

    Boundary(x_start, x_end)

The boundary is a continuous line segment `Segment` between `x_start` and `x_end`.

### 2D Constructors

    Boundary(verts)

The boundary is a closed polygon `Ngon` with vertices `verts`, which should be be a
vector of `(x, y)` coordinates or a `2 x N` matrix.

    Boundary(xmin, xmax, ymin, ymax)

A convenience constructor is also provided for the case of a rectangular boundary.

### ≥3D Constructors

    Boundary(corner_min, corner_max)

The boundary is a `HyperRectangle` with minimum and maximum vertices `corner_min`, `corner_max`
where the corners are `NTuple`s, `(x_min, y_min, z_min, ...)`, `(x_max, y_max, z_max, ...)`.
"""
struct Boundary{Dim, M, CRS}
    boundary::Polytope{Dim, M, CRS}

    function Boundary(; boundary::Polytope{Dim, M, CRS}) where {Dim, M, CRS}
        return new{Dim, M, CRS}(boundary)
    end
end

### 1D
function Boundary(x_start::Real, x_end::Real)
    @argcheck x_start < x_end

    return Boundary(boundary = Segment(Meshes.Point(float(x_start)), Meshes.Point(float(x_end))))
end

function points(boundary::Boundary{1, M, CRS}) where {M, CRS}
    verts = boundary.boundary.vertices
    
    return stack([[coords(v).x.val] for v in verts])
end

### 2D
function Boundary(verts::Vector{<:Tuple{Real, Real}})
    @argcheck all(length.(verts) .== 2)

    return Boundary(boundary = Ngon(verts...))
end

function Boundary(verts::Matrix{<:Real})
    @argcheck size(verts, 1) == 2

    return Boundary([(verts[1, i], verts[2, i]) for i = 1:size(verts, 2)])
end

function Boundary(xmin::Real, xmax::Real, ymin::Real, ymax::Real)
    @argcheck xmin < xmax
    @argcheck ymin < ymax

    return Boundary([(xmin, ymin), (xmax, ymin), (xmax, ymax), (xmin, ymax)])
end

function points(boundary::Boundary{2, M, CRS}) where {M, CRS}
    verts = boundary.boundary.vertices
    
    return stack([[coords(v).x.val, coords(v).y.val] for v in verts])
end

### ≥3D
function Boundary(corner_min::NTuple{N, Real}, corner_max::NTuple{N, Real}) where {N}
    @argcheck length(corner_min) == length(corner_max)
    @argcheck N ≥ 3 "Prefer to not use a HyperRectangle in 1 or 2 dimensions"
    return Boundary(boundary = HyperRectangle(corner_min, corner_max))
end

### GENERAL

# autoboundary
## start with bounding box of data
## iteratively shrink (bisect ?) until a certain fraction of data is in nirvana