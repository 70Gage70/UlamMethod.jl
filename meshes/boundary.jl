"""
    struct Boundary{K, Dim, CRS}

A container type for the computational domain within which Ulam's method is applied.

The boundary is partitioned according to an [`BinningAlgorithm`](@ref).

### Fields

- `boundary`: A `Polytope` defining the boundary.

### 1D Constructor

In 1D, the boundary is a continuous line segment between `x_start` and `x_end`. Use

`Boundary(x_start, x_end)`

### 2D Constructors

In 2D, the boundary is a closed polygon with vertices `verts`, which should be be a
vector of `[x, y]` coordinates. Use

`Boundary(verts)`

A convenience constructor is also provided for the case of a rectangular boundary. Use

`Boundary(xmin, xmax, ymin, ymax)`
"""
struct Boundary{K, Dim, CRS}
    boundary::Polytope{K, Dim, CRS}

    function Boundary(; boundary::Polytope{K, Dim, CRS} = boundary) where {K, Dim, CRS}
        return new{K, Dim, CRS}(boundary)
    end
end

### 1D
function Boundary(x_start::Real, x_end::Real)
    @argcheck x_start < x_end

    return Boundary(boundary = Segment(Meshes.Point(float(x_start)), Meshes.Point(float(x_end))))
end

### 2D
function Boundary(verts::Vector{<:Tuple{Real, Real}})
    @argcheck all(length.(verts) .== 2)

    return Boundary(boundary = Ngon(verts...))
end

function Boundary(xmin::Real, xmax::Real, ymin::Real, ymax::Real)
    @argcheck xmin < xmax
    @argcheck ymin < ymax

    return Boundary([(xmin, ymin), (xmax, ymin), (xmax, ymax), (xmin, ymax)])
end