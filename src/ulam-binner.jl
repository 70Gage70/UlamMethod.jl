
include("binner-voronoi.jl")
include("binner-hexagon.jl")
include("binner-square.jl")

"""
    ulam_binner(traj, domain)

Selects and executes the appropriate binning algorithm for the given trajectories and domain.
Intersects the resulting polygons with the domain so that the returned polygons are clipped to it.
"""
function ulam_binner(traj::UlamTrajectories, domain::UlamDomain)::Vector{UlamPolygon}
    poly_type = domain.poly_type

    if poly_type == "reg"
        res = binner_square(domain)
    elseif poly_type == "hex"
        res = binner_hexagon(domain)
    elseif poly_type == "vor"        
        res = binner_voronoi(traj, domain)
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