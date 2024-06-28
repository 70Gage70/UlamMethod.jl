# Meshes Formulation

- The user should provide a shape that defines the boundary.
    - 1D: `<:Chain`, i.e. 1D `Polytope`
    - 2D: `<:Polygon`, i.e. 2D `Polytope`

- `Bins` is a struct with one field `bins` that holds the bins; it is a vector of `Geometry`s, essentially equivalent to a `GeometrySet`.

- `AbstractBinner{Dim}` is the abstract type for a binning algorithm in dimension `Dim`. Each subtype should implement a method `_bin(boundary, binner)` which returns a `Bins`
    - 1D: `LineBinner <: AbstractBinner{1}`
    - 1D: `RectangleBinner <: AbstractBinner{2}`

- `bin` handles the dispatch to the appropriate `_bin` and any other high level binning options.

- `Trajectories{Dim}` holds `x0` and `xT` as `Dim x n_points` matrices.

