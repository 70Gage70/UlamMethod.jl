
include("binner-voronoi.jl")
include("binner-hexbin.jl")
include("binner-square.jl")

"""
    All binners.
"""

function ulam_binner(traj::UlamTrajectories, domain::UlamDomain)::UlamCovering
    bin_type = domain.bin_type

    if bin_type == "reg"
        res = square_binner(traj, domain) # rseed = 123 is "default"
    elseif bin_type == "hex"
        # res = hexbin_binner(n_polys, corners) 
        error("Not ready yet.") 
    elseif bin_type == "vor"        
        # res = voronoi_binner(n_polys, corners)
        error("Not ready yet.") 
    end

    return res
end