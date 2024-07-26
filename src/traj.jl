"""
    struct Trajectories{Dim}

A container for trajectory data of dimension `Dim`.

### Fields

- `x0`: Initial coordinates, should be a matrix of size `Dim` x `n_traj`.
- `xT`: Final coordinates, should be a matrix of size `Dim` x `n_traj`.

### Constructors

    Trajectories(x0, xT)

Construct `Trajectories` directly from data.

    Trajectories(dim, n_traj; corners = nothing)

Construct `n_traj` random `Trajectories` of dimension `dim` for testing. The trajectories will lie inside the box defined 
by `corners = ((xmin, ymin, zmin, ...), (xmax, ymax, zmax, ...))` which defaults to the unit box when not provided.
"""
struct Trajectories{Dim}
    x0::Matrix{Float64}
    xT::Matrix{Float64}

    function Trajectories(x0::Matrix{<:Real}, xT::Matrix{<:Real})
        @argcheck size(x0) == size(xT)
        @argcheck all(size(x0) .> 0)

        return new{size(x0, 1)}(x0, xT)
    end
end

function Trajectories(
    dim::Integer, 
    n_traj::Integer; 
    corners::Union{Tuple{NTuple{N, Real}, NTuple{N, Real}}, Nothing} = nothing) where {N}

    if corners === nothing
        corners = (zeros(dim), ones(dim))
    else
        @argcheck dim == N
    end

    box_min, box_max = corners
    @argcheck all(box_max .> box_min)

    x0 = [box_min[i] + (box_max[i] - box_min[i])*rand() for i = 1:dim, _ = 1:n_traj]
    xT = [box_min[i] + (box_max[i] - box_min[i])*rand() for i = 1:dim, _ = 1:n_traj]

    return Trajectories(x0, xT)
end