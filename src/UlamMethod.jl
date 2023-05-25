module UlamMethod

export ulam_method, ulam_write # from main.jl
export # from types.jl
    UlamPolygon,
    UlamTrajectories,
    UlamDomain,
    UlamInfo,
    UlamResult,
    PolyTable,
    InpolyResult,
    show
export ulam_intersection, ulam_intersects, inpoly, ulam_polys_to_indices # from helpers.jl


include("types.jl")
include("main.jl")

module EarthUlamPolygons
    include("earth-polygons.jl")
end

end # module