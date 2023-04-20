import MAT
import HDF5
import LinearAlgebra

"""
    global_poly_types

One of `"sqr"`, `"hex"` and `"vor"` according to the coverings in [`UlamDomain`](@ref).
"""
const global_poly_types::Vector{String} = ["sqr", "hex", "vor"]

"""
    global_stoc_types

One of `"data"` or `"source"` according to the stochasticization options in [`UlamDomain`](@ref).
"""
const global_stoc_types::Vector{String} = ["data", "source"]

"""
    global_traj_file_types

One of `"mat"`, `"h5"` according to the trajectory input file in [`UlamTrajectories`](@ref)
"""
const global_traj_file_types::Vector{String} = ["mat", "h5"]

"""
    global_poly_number_default

The default number of polygons for each of the types of coverings.

- `"sqr"`: 500
- `"hex"`: 500
- `"vor"`: 100
"""
const global_poly_number_default::Dict = Dict("sqr" => 500, "hex" => 500, "vor" => 100)

"""
    global_rseed_default

The default seed for RNG reproducibility. Default `123`.
"""
const global_rseed_default::Int64 = 123

include("types/UlamPolygon.jl")
include("types/UlamTrajectories.jl")
include("types/UlamDomain.jl")
include("types/UlamInfo.jl")
include("types/UlamResult.jl")
include("types/PolyTable.jl")
include("types/InpolyResult.jl")