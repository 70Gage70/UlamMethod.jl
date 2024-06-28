"""
    abstract type ReinjectionAlgorithm

The abstract type for the algorithm controlling how trajectories pointing from nirvana to the interior are redistributed.

Each subtype `reinj_algo` should implement the following:

- a function `_reinject(bins, Pij, reinj_algo)` that computes a new transition matrix `Pij` given the old one 
and the `bins`. This function is applied after the largest strongly connected component of the original `Pij`
is taken, but before the matrix is row-stochasticized. 
"""
abstract type ReinjectionAlgorithm end

"""
    struct DataReinjection

Reinject the data according to which bins are actually hit by transitions from the exterior to the interior. 

Since this is applied by default in the main Ulam's method algorithm, the associated `_reinject` function 
for `DataReinjection` just returns the input `Pij`.
"""
struct DataReinjection <: ReinjectionAlgorithm
end

function _reinject(bins::Bins, Pij::Matrix{Float64}, reinj_algo::DataReinjection)
    return Pij
end

"""
    reinject(bins, Pij, reinj_algo)

Compute the new transition matrix `Pij` according to the old one and `reinj_algo::ReinjectionAlgorithm`, and 
the bins that contain data.
"""
function reinject(bins::Bins, Pij::Matrix{Float64}, reinj_algo::ReinjectionAlgorithm)
    return _reinject(bins, Pij, reinj_algo)
end