"""
    struct Trajectories{Dim}

A container for trajectory data of dimension `Dim`.

### Fields

- `x0`: Initial coordinates, should be a matrix of size `Dim` x `n_traj`.
- `xT`: Final coordinates, should be a matrix of size `Dim` x `n_traj`.

### Constructors

    Trajectories(x0, xT)

Construct `Trajectories` directly from data.

    Trajectories(dim, n_traj; corners = nothing, mu_sigma = nothing)

Construct `n_traj` random `Trajectories` of dimension `dim` for testing. 

The `x0` will lie inside the box defined by `corners = ((xmin, ymin, zmin, ...), (xmax, ymax, zmax, ...))` which 
defaults to the unit box when not provided.

The `xT` are displaced from the `x0` by normal distributions along each dimension with means and variances defined by 
`mu_sigma = [(mu1, sigma1), (mu2, sigma2), ...]`. If not provided, default to `mu = 1.0`, `sigma = 1.0`
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
    corners::Union{Tuple{NTuple{N, Real}, NTuple{N, Real}}, Nothing} = nothing,
    mu_sigma::Union{Vector{<:Tuple{Real, Real}}, Nothing} = nothing) where {N}

    if corners === nothing
        corners = (zeros(dim), ones(dim))
    else
        @argcheck dim == N
    end

    if mu_sigma === nothing
        mu_sigma = [(1.0, 1.0) for _ = 1:dim]
    else
        @argcheck dim == length(mu_sigma)
    end

    box_min, box_max = corners
    @argcheck all(box_max .> box_min)

    x0 = stack([rand(Distributions.Uniform(box_min[i], box_max[i]), n_traj) for i = 1:dim], dims = 1)
    xT = x0 + stack([rand(Distributions.Normal(mu_sigma[i]...), n_traj) for i = 1:dim], dims = 1)

    return Trajectories(x0, xT)
end

Base.length(traj::Trajectories) = size(traj.x0, 2)