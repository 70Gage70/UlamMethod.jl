"""
General helper functions.
"""
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

using .UlamTypes

import MAT
import HDF5
import PolygonInbounds

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################



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
    # findall(x->x==true, ip2res[:,1,i]) therfore finds the indices of all data points in polygon i
    # contains is a lookup table for polygons that contain data, 
    # e.g. contains[6] = 3 means that the 6th polygon is the 3rd polygon in the list that contains data
    # note that [i for i in keys(contains)] is a list of indices of polygons that contain data

    return InpolyResult(inds, contains)
end

"""
    ulam_intersect(verts1, verts2)

Compute the intersection of two `UlamPolygons` objects.
"""

function ulam_intersect(verts1::Matrix{<:Real}, verts2::Matrix{<:Real})
    if typeof(verts1) != Matrix{Float64}
        verts1 = convert(Matrix{Float64}, verts1)
    end

    if typeof(verts2) != Matrix{Float64}
        verts2 = convert(Matrix{Float64}, verts2)
    end

    if verts1[end, :] != verts1[1, :]
        verts1 = vcat(verts1, verts1[1,:]')
    end

    if verts2[end, :] != verts2[1, :]
        verts2 = vcat(verts2, verts2[1,:]')
    end    

    verts2 = convert(Matrix{Float64}, verts2)
    p1 = LibGEOS.Polygon([[verts1[i,:] for i = 1:size(verts1, 1)]])
    p2 = LibGEOS.Polygon([[verts2[i,:] for i = 1:size(verts2, 1)]])
    pint = GeoInterface.coordinates(LibGEOS.intersection(p1, p2))

    if length(pint[1]) == 0
        return false
    end

    res = Vector{Matrix{Float64}}()
    
    for p in pint
        push!(res, reduce(hcat, p[1])')
    end

    return res    
end



function help_smear(xmin, xmax, ymin, ymax; resolution = 0.1)
    centers = collect.(Iterators.product(xmin:resolution:xmax, ymin:resolution:ymax))
    centers = reshape(centers, (length(centers),1))
    centers = transpose(reduce(hcat, centers))
    return centers
end


function ABcorner_standard()
    corners = [-100, 15, -9, 39]

    A_centers = [
        -18.0 17.0;
        ]
    B_centers = help_smear(-98.0, -92.0, 17.7, 32.0, resolution = 0.1)

    return corners, A_centers, B_centers
end

function corners_mid(corners)
    return [(corners[1] + corners[2])/2, (corners[3] + corners[4])/2]
end

# takes in a vector of vectors and converts it to a matrix
function vecvec_to_mat(vvec)
    return reduce(hcat,vvec)'
end

function write_dict_to_h5(fout, group_name, dict)
    g = create_group(fout, group_name)
    for key in collect(keys(dict))
        if dict[key] == []
            # g[key] = HDF5.EmptyArray{Float64}()
            g[key] = zeros(0)
        else
            g[key] = dict[key]
        end
    end
end
