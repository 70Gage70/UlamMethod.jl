
include("binners/binner-voronoi.jl")
include("binners/binner-hexagon.jl")
include("binners/binner-square.jl")
include("binners/binner-rectangle.jl")
include("binners/binner-triangle.jl")

"""
    ulam_binner(traj, domain)

Select and executes the appropriate binning algorithm for the given trajectories and domain.

Intersect the resulting polygons with the domain so that the returned polygons are clipped to it.
"""
function ulam_binner(traj::UlamTrajectories, domain::UlamDomain)
    poly_type = domain.poly_type

    if poly_type == "sqr"
        res = binner_square(domain)
    elseif poly_type == "rec"
        res = binner_rectangle(domain)        
    elseif poly_type == "hex"
        res = binner_hexagon(domain)
    elseif poly_type == "tri"
        res = binner_triangle(domain)
    elseif poly_type == "vor"        
        res = binner_voronoi(traj, domain)
    end

    # intersect the resulting polygons with the domain
    # the polygons are now exactly clipped to the boundary of nirvana by the time they get to ulam_nirvana
    intersected = Vector{UlamPolygon{Float64}}()
    for poly in res
        int = ulam_intersection(domain.domain, poly)
        if int != false
            push!(intersected, int)
        end
    end   
    
    return intersected
end