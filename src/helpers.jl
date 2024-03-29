import PolygonInbounds
import LibGEOS, GeoInterface

"""
    inpoly(data::Matrix{<:Real}, polys::PolyTable)

Determine which polygon of `polys` contains the data points in `data`. Return an [`InpolyResult`](@ref).
"""
function inpoly(data::Matrix{<:Real}, polys::PolyTable)
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

Determine the indices of points that are inside `domain.domain`. Return an [`InpolyResult`](@ref).
"""
function inpoly(traj::UlamTrajectories, domain::UlamDomain)
    data = [traj.x0 ;; traj.y0]
    polys = PolyTable([domain.domain])

    return inpoly(data, polys)
end


"""
    ulam_intersection(poly1, poly2)

Compute the intersection of two [`UlamPolygon`](@ref) objects. Return `false` if they do not intersect or if they only intersect along a point or line.
"""
function ulam_intersection(poly1::UlamPolygon, poly2::UlamPolygon)
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
    pkind = GeoInterface.geomtrait(pint)

    # Ensure that the result is a polygon
    if pkind == GeoInterface.PolygonTrait()
        pint = GeoInterface.coordinates(pint)

        if length(pint[1]) == 0
            return false # no intersection
        else
            pmat = reduce(hcat, pint[1])'
            pmat = pmat[1:end-1,:]

            return UlamPolygon(pmat)
        end
    elseif pkind == GeoInterface.MultiPolygonTrait()
        # there are multiple polygons in the intersection; take the largest one
        all_ints = GeoInterface.getgeom(pint)
        largest_p = argmax([LibGEOS.area(p) for p in all_ints])
        pint = GeoInterface.coordinates(collect(all_ints)[largest_p])

        pmat = reduce(hcat, pint[1])'
        pmat = pmat[1:end-1,:]

        return UlamPolygon(pmat)
    elseif pkind in [GeoInterface.LineStringTrait(), GeoInterface.PointTrait()]      
        return false # polygons only intersect at a point, or along a line  
    else
        @error "The intersection of these polygons is $(pkind), this shouldn't happen." poly1.nodes poly2.nodes 
    end
end

"""
    ulam_intersects(poly1, poly2)

Return `true` or `false` accoring to whether two [`UlamPolygon`](@ref) objects intersect.

This is a call to [`ulam_intersection`](@ref), hence two polygons are not considered to intersect if their intersection is only along a point or line.
"""
function ulam_intersects(poly1::UlamPolygon, poly2::UlamPolygon)
    res = ulam_intersection(poly1, poly2)

    if res == false
        return false
    else
        return true
    end
end

"""
    ulam_polys_to_indices(ulam_res::Vector{<:UlamPolygon}, region)

Find the indicies of `polys` which contain `region`.

- If `region` is entered as an `N x 2` matrix of numbers, then the indices are of polygons which contain at least one point from `region`.
- If `region` is entered as an `UlamPolygon`, then the indices are of polygons which intersect with `region`.
"""
function ulam_polys_to_indices(polys::Vector{<:UlamPolygon}, region::Union{Matrix{<:Real}, UlamPolygon})
    if region isa UlamPolygon
        inds = [i for i in 1:length(polys) if ulam_intersects(region, polys[i])]
    elseif region isa Matrix
        inds = unique(inpoly(region, PolyTable(polys)).inds)
    end

    filter!(x->x!=0, inds)

    return inds
end

"""
    ulam_polys_to_indices(ulam_res::UlamResult, region)

When applied to `UlamResult`, find the container indices of `ulam_res.polys`.
"""
function ulam_polys_to_indices(ulam_res::UlamResult, region::Union{Matrix{<:Real}, UlamPolygon})
    polys = ulam_res.polys

    if region isa UlamPolygon
        inds = [i for i in 1:length(polys) if ulam_intersects(region, polys[i])]
    elseif region isa Matrix
        inds = unique(inpoly(region, PolyTable(polys)).inds)
    end

    filter!(x->x!=0, inds)

    return inds
end

