"""
General helper functions.
"""
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

using .UlamTypes

import MAT, HDF5
import PolygonInbounds
import LibGEOS, GeoInterface

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################


"""
    inpoly(data::Matrix{Float64}, polys::PolyTable)

Determines which polygon of `polys` contains the data points in `data`. Returns an `InpolyResult`[@ref].
"""
function inpoly(data::Matrix{Float64}, polys::PolyTable)::InpolyResult
    @assert size(data, 1) > 0
    @assert size(data, 2) == 2

    ip2res = PolygonInbounds.inpoly2(data, polys.nodes, polys.edges)
    inds = zeros(Int64, size(data, 1))
    contains = Dict{Int64, Int64}()

    j = 1
    for i = 1:size(ip2res, 3)
        fa = findall(x->x==true, ip2res[:,1,i])
        if length(fa) > 0
            contains[i] = j
            j = j + 1
            inds[fa] .= i
        end
    end

    # ip2res is a BitArray with dimensions size(data, 1) x 2 x size(polys.nodes, 1).
    # ip2res[:,1,i] is a BitVector such that ip2res[:,1,i][k] == true if the k'th data point is in polygon i
    # fa = findall(x->x==true, ip2res[:,1,i]) therefore finds the indices of all data points in polygon i
    # contains is a lookup table for polygons that contain data, 
    # e.g. contains[6] = 3 means that the 6th polygon is the 3rd polygon in the list that contains data
    # note that [i for i in keys(contains)] is a list of indices of polygons that contain data

    return InpolyResult(inds, contains)
end

"""
    inpoly(traj::UlamTrajectories, domain::UlamDomain)

Determines the indices of points that are inside the domain (i.e. not nirvana.)
"""
function inpoly(traj::UlamTrajectories, domain::UlamDomain)::InpolyResult
    data = [traj.x0 ;; traj.y0]

    if domain.domain === nothing # the domain is a square
        xmin, xmax, ymin, ymax = domain.corners
        polys = PolyTable([UlamPolygon([xmin ymin; xmin ymax; xmax ymax; xmax ymin])])
    else
        polys = PolyTable([domain.domain])
    end

    return inpoly(data, polys)
end


"""
    ulam_intersection(poly1, poly2)

Compute the intersection of two `UlamPolygons` objects. Returns false if they do not intersect.
"""
function ulam_intersection(poly1::UlamPolygon, poly2::UlamPolygon)::Union{Bool,UlamPolygon}
    # UlamPolygons are not closed, so we have to close them for LibGEOS.
    nodes1 = poly1.nodes
    nodes2 = poly2.nodes
    nodes1 = vcat(nodes1, nodes1[1,:]')
    nodes2 = vcat(nodes2, nodes2[1,:]')

    # Construct LibGEOS-formatted polygons
    p1 = LibGEOS.Polygon([[nodes1[i,:] for i = 1:size(nodes1, 1)]])
    p2 = LibGEOS.Polygon([[nodes2[i,:] for i = 1:size(nodes2, 1)]])

    # Compute the intersection
    pint = LibGEOS.intersection(p1, p2)

    # Ensure that the result is a polygon
    if GeoInterface.geomtrait(pint) == GeoInterface.PolygonTrait()
        pint = GeoInterface.coordinates(pint)

        if length(pint[1]) == 0
            return false # no intersection
        else
            pmat = reduce(hcat, pint[1])'
            pmat = pmat[1:end-1,:]

            return UlamPolygon(pmat)
        end
    elseif GeoInterface.geomtrait(pint) == GeoInterface.MultiPolygonTrait()
        @error "The intersection of these polygons is a multipolygon, this shouldn't happen." poly1 poly2
    else
        @error "The intersection of these polygons is bad, this REALLY shouldn't happen." poly1 poly2
    end
end

"""
    ulam_intersects(poly1, poly2)

Return truw or false accoring to whether two `UlamPolygons` objects intersect. Can be faster
than `ulam_intersection` if the shape of the intersection is not needed.
"""
function ulam_intersects(poly1::UlamPolygon, poly2::UlamPolygon)::Bool
    # UlamPolygons are not closed, so we have to close them for LibGEOS.
    nodes1 = poly1.nodes
    nodes2 = poly2.nodes
    nodes1 = vcat(nodes1, nodes1[1,:]')
    nodes2 = vcat(nodes2, nodes2[1,:]')

    # Construct LibGEOS-formatted polygons
    p1 = LibGEOS.Polygon([[nodes1[i,:] for i = 1:size(nodes1, 1)]])
    p2 = LibGEOS.Polygon([[nodes2[i,:] for i = 1:size(nodes2, 1)]])

    # Compute the intersection
    return LibGEOS.intersects(p1, p2)
end