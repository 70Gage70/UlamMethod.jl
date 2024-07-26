"""
    struct Bins{Dim, M, CRS}

A container type for bins of dimension `Dim` embedded in manifold `M` with coordinate reference system `CRS`.

### Fields

- `bins`: A vector of `Polytope{Dim, M, CRS}` objects.

### Methods

    points(bins)
    
Return a vector of raw vertices for each bin.
"""
struct Bins{Dim, M, CRS}
    bins::Vector{<:Polytope{Dim, M, CRS}}
end

"""
    points(bins)

Return a vector of raw vertices for each bin. Can also be applied to a [`BinningAlgorithm`](@ref).
"""
function points(bins::Bins{Dim, M, CRS}) where {Dim, M, CRS}
    bins = bins.bins

    if Dim == 1
        return [[(c.x.val) for c in coords.(bin.vertices)] for bin in bins]
    elseif Dim == 2
        return [[(c.x.val, c.y.val) for c in coords.(bin.vertices)] for bin in bins]
    else
        [vec(collect(Iterators.product([(bin.min[i], bin.max[i]) for i = 1:length(bin.min)]...))) for bin in bins]
    end
end

"""
    abstract type BinningAlgorithm{Dim}

An abstract type for binning algorithms of dimension `Dim`.

### Implementation

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

### Methods

    points(binner)
    
Return a vector of raw vertices for each bin.
"""
abstract type BinningAlgorithm{Dim} end

points(binner::BinningAlgorithm) = points(binner.bins)