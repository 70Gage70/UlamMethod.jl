"""
    struct Boundary{Dim, CRS}

A container type for the computational domain within which Ulam's method is applied.

The boundary is partitioned according to an [`BinningAlgorithm`](@ref).

### Fields

- `boundary`: A `Polytope` defining the boundary.

### 1D

The boundary is a continuous line segment `Segment` between `x_start` and `x_end`. Use

`Boundary(x_start, x_end)`

Use `points(boundary)` to calculate a `Dim x N` matrix of boundary vertices.

### 2D 

The boundary is a closed polygon `Ngon` with vertices `verts`, which should be be a
vector of `(x, y)` coordinates or a `2 x N` matrix. Use

`Boundary(verts)`

A convenience constructor is also provided for the case of a rectangular boundary. Use

`Boundary(xmin, xmax, ymin, ymax)`

Use `points(boundary)` to calculate a `Dim x N` matrix of boundary vertices.

### ≥3D

The boundary is a `HyperRectangle` with minimum and maximum vertices `corner_min`, `corner_max`. Use

`Boundary(corner_min, corner_max)`

where the corners are `NTuple`s, `(x_min, y_min, z_min, ...)`, `(x_max, y_max, z_max, ...)`.
"""
struct Boundary{Dim, CRS}
    boundary::Polytope{Dim, Dim, CRS}

    function Boundary(; boundary::Polytope{Dim, Dim, CRS}) where {Dim, CRS}
        return new{Dim, CRS}(boundary)
    end
end

### 1D
function Boundary(x_start::Real, x_end::Real)
    @argcheck x_start < x_end

    return Boundary(boundary = Segment(Meshes.Point(float(x_start)), Meshes.Point(float(x_end))))
end

function points(boundary::Boundary{1, CRS}) where {CRS}
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

function points(boundary::Boundary{2, CRS}) where {CRS}
    verts = boundary.boundary.vertices
    
    return stack([[coords(v).x.val, coords(v).y.val] for v in verts])
end

### ≥3D
function Boundary(corner_min::NTuple{N, Real}, corner_max::NTuple{N, Real}) where {N}
    @argcheck length(corner_min) == length(corner_max)
    @argcheck N ≥ 3 "Prefer to not use a HyperRectangle in 1 or 2 dimensions"
    return Boundary(boundary = HyperRectangle(corner_min, corner_max))
end