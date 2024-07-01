"""
    struct HyperRectangle{K, Dim, CRS}

An extention of `Meshes` to dimensions greater than 3.

The extension is minimal and implements no methods and does not interface with CRS.

### Fields 

- `min`: The coordinates of the minimum of the HyperRectangle as an `NTuple`.
- `max`: The coordinates of the maximum of the HyperRectangle as an `NTuple`.
- `vertices`: `[min, max]`, required by `Meshes`.

### Constructor

`HyperRectangle(min, max)`
"""
struct HyperRectangle{K, Dim, CRS} <: Polytope{K, Dim, CRS}
    min::NTuple{Dim, Float64}
    max::NTuple{Dim, Float64}
    vertices::Vector{NTuple{Dim, Float64}}

    function HyperRectangle(min::NTuple{Dim, Real}, max::NTuple{Dim, Real}) where {Dim}
        t = typeof(Box(tuple(zeros(Dim)...), tuple(zeros(Dim)...))).parameters[2]
        return new{Dim, Dim, t}(min, max, [min, max])
    end
end