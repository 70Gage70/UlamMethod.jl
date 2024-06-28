"""
    struct Boundary{K, Dim, CRS}

A container type for the computational domain within which Ulam's method is applied.

The boundary is partitioned according to an [`AbstractBinner`](@ref).

### Fields

- `boundary`: A `Polytope` defining the boundary.

### Constructor

- 1D: `Boundary(boundary::Segment)`
- 2D: `Boundary(boundary::Ngon)`
"""
struct Boundary{K, Dim, CRS}
    boundary::Polytope{K, Dim, CRS}

    function Boundary(; boundary::Polytope{K, Dim, CRS} = boundary) where {K, Dim, CRS}
        return new{K, Dim, CRS}(boundary)
    end
end

### 1D
Boundary(boundary::Segment) = Boundary(boundary = boundary)

### 2D
Boundary(boundary::Ngon) = Boundary(boundary = boundary)