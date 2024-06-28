"""
    struct Boundary{K, Dim, CRS}

A container type for the computational domain within which Ulam's method is applied.

The boundary is partitioned according to an [`AbstractBinner`](@ref).

Data outside the boundary are considered in nirvana regardless of binning.

### Fields

- `boundary`: A `Polytope` whose dimension matches the dimension of the data.

### Constructor

`Boundary(boundary)`
"""
struct Boundary{Dim, CRS}
    boundary::Polytope{Dim, Dim, CRS}

    function Boundary(boundary::Polytope{Dim, Dim, CRS}) where {Dim, CRS}
        new{Dim, CRS}(boundary)
    end
end