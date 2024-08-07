"""
    abstract type ReinjectionAlgorithm

The abstract type for the algorithm controlling how trajectories pointing from nirvana to the interior are redistributed.

Each subtype `reinj_algo` should implement the following:

- a function `reinject!(binner, Pij, reinj_algo)` that modifies the transition matrix `Pij`in place given
the binning algorithm `binner`. This function is applied after the largest strongly connected component of the original `Pij`
is taken, but before the matrix is row-stochasticized. 
"""
abstract type ReinjectionAlgorithm end

"""
    struct DataReinjection

Reinject the data according to which bins are actually hit by transitions from the exterior to the interior. 

Since this is applied by default in the main Ulam's method algorithm, the associated `reinject!` function 
for `DataReinjection` does nothing.

### Constructor

    DataReinjection()
"""
struct DataReinjection <: ReinjectionAlgorithm
end

function reinject!(binner::BinningAlgorithm, Pij::Matrix{Float64}, reinj_algo::DataReinjection)
    return nothing
end

"""
    struct SourceReinjection{Dim}

Reinject the data at particular locations. 

### Fields

- `points`: A vector of the form `[p1, p2, ...]` where `pi` is a point such as `(x)` or `(x, y)`. Data are \
reinjected uniformly across all bins that contain at least one member of `points`.
- `fallback`: The `ReinjectionAlgorithm` to apply in case no bins contain any points or there is no data to reinject.

### Constructor

    SourceReinjection(points; fallback = DataReinjection())
"""
struct SourceReinjection{Dim} <: ReinjectionAlgorithm
    points::Vector{NTuple{Dim, Float64}}
    fallback::ReinjectionAlgorithm

    function SourceReinjection(
        points::Vector{<:NTuple{N, Real}}; 
        fallback::ReinjectionAlgorithm = DataReinjection()) where {N}

        return new{N}(points, fallback)
    end
end

function reinject!(
    binner::BinningAlgorithm{Dim}, 
    Pij::Matrix{Float64}, 
    reinj_algo::SourceReinjection{Dim}) where {Dim}

    memb = membership(stack(collect.(reinj_algo.points)), binner)
    total_counts = sum(Pij[end,:])

    if all(isnothing.(memb)) || total_counts == 0
        @warn "Could not reinject at source - using fallback."
        return reinject!(binner, Pij, reinj_algo.fallback)
    end

    memb = memb[findall(!isnothing, memb)]
    Pij[end, :] .= 0
    Pij[end, memb] .= 1
    return nothing
end

"""
    struct StationaryReinjection

Reinject the data weighted by the stationary distribution of `P_open`.

### Fields

- `fallback`: The `ReinjectionAlgorithm` to apply in case the stationary distribution can not be computed.

### Constructor

    StationaryReinjection(; fallback = DataReinjection())
"""
struct StationaryReinjection <: ReinjectionAlgorithm
    fallback::ReinjectionAlgorithm

    function StationaryReinjection(; fallback::ReinjectionAlgorithm = DataReinjection())
        return new(fallback)
    end
end

function reinject!(
    binner::BinningAlgorithm{Dim}, 
    Pij::Matrix{Float64}, 
    reinj_algo::StationaryReinjection) where {Dim}

    try
        p = Pij[1:end-1, 1:end-1]
        Pij[end,1:end-1] .= abs.(eigvecs(p')[:,end])
    catch
        @warn "Could not reinject using stationary distribution - using fallback."
        return reinject!(binner, Pij, reinj_algo.fallback)
    end

    return nothing
end