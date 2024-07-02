"""
    struct Bins{Dim, CRS}

A container type for bins.

### Fields

- `bins`: A vector of `Polytope` objects of dimension `Dim` and coordinate reference system `CRS`.
"""
struct Bins{Dim, CRS}
    bins::Vector{<:Polytope{Dim, Dim, CRS}}
end

"""
    abstract type BinningAlgorithm{Dim}

An abstract type for binning algorithms of dimension `Dim`.

Each subtype should have a fields:

- `boundary::Boundary{Dim, CRS}``
- `bins::Bins{Dim, CRS}`
- `idx2pos::Vector{Union{Int64, Nothing}}`

Each subtype should implement the following:

- A function `membership(data::Matrix, binner)` which takes a `Dim x N_points` matrix of points and \
returns a vector `memb` where `memb[i] = j` if `data[:,i]` is inside `binner.bins[j]` and `memb[i] = nothing` \
if `data[:,i]` is not inside any bin.
- A function `membership(traj::Trajectories, binner)` which returns `(memb_x0, memb_xT)`. In many cases this \
will be given by `(membership(traj.x0, binner), membership(traj.xT, binner))`.
"""
abstract type BinningAlgorithm{Dim} end