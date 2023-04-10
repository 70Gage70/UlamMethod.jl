module UlamMethod

export ulam_method # from main.jl
export # from types.jl
    UlamPolygon,
    UlamTrajectories,
    UlamDomain,
    UlamInfo,
    UlamResult,
    PolyTable,
    InpolyResult,
    show

include("types.jl")
include("main.jl")

end # module