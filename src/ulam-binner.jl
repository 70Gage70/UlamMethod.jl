
include("binner-voronoi.jl")
include("binner-hexbin.jl")
include("binner-square.jl")

"""
    All binners.
"""

function ulam_binner(traj::UlamTrajectories, domain::UlamDomain)::Vector{UlamPolygon}
    bin_type = domain.bin_type

    if bin_type == "reg"
        res = square_binner(traj, domain)
    elseif bin_type == "hex"
        # res = hexbin_binner(n_polys, corners) 
        error("Not ready yet.") 
    elseif bin_type == "vor"        
        # res = voronoi_binner(n_polys, corners)
        error("Not ready yet.") 
    end

    # intersect the binned polygons with the given domain if one exists
    if domain.domain != nothing
        intersected = Vector{UlamPolygon}()
        for poly in res
            int = ulam_intersection(domain, poly)
            if int != false
                push!(intersected, int)
            end
        end

        res = intersected
    end
    
    return res
end