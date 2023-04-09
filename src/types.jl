module UlamTypes

import MAT
import HDF5
import LinearAlgebra

export 
    UlamPolygon,
    UlamTrajectories,
    UlamDomain,
    UlamInfo,
    UlamResult,
    PolyTable,
    InpolyResult,
    show

const global_poly_types::Vector{String} = ["sqr", "hex", "vor"]
const global_stoc_types::Vector{String} = ["data", "source"]
const global_traj_file_types::Vector{String} = ["mat", "h5"]
const global_poly_number_default::Dict = Dict("sqr" => 500, "hex" => 500, "vor" => 100)

include("types/UlamPolygon.jl")
include("types/UlamTrajectories.jl")
include("types/UlamDomain.jl")
include("types/UlamInfo.jl")
include("types/UlamResult.jl")
include("types/PolyTable.jl")
include("types/InpolyResult.jl")

end # module

