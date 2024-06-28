"""
    struct Trajectories{Dim}

A container for trajectory data of dimension `Dim`.

### Fields

- `x0`: Initial coordinates, should be a matrix of size `Dim` x `n_traj`.
- `xT`: Final coordinates, should be a matrix of size `Dim` x `n_traj`.

### Constructor

`Trajectories(x0, xT)`
"""
struct Trajectories{Dim}
    x0::Matrix{Float64}
    xT::Matrix{Float64}

    function Trajectories(x0::Matrix{<:Real}, xT::Matrix{<:Real})
        @argcheck size(x0) == size(xT)
        @argcheck all(size(x0) .> 0)
        @argcheck size(x0, 1) in [1, 2] "Dimension of data must be 1 or 2."

        return new{size(x0, 1)}(x0, xT)
    end
end
