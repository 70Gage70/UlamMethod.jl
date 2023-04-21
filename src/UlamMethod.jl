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
export ulam_intersection, ulam_intersects, inpoly # from helpers.jl
export North_Atlantic_clipped_verts, North_Atlantic_box_verts, GoG_big_verts, GoG_small_verts # from earth_polygons.jl

include("types.jl")
include("main.jl")
include("earth-polygons.jl")

end # module