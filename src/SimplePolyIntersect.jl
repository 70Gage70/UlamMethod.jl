module SimplePolyIntersect

import LibGEOS
import GeoInterface

export polyintersect

function polyintersect(verts1::Matrix{<:Real}, verts2::Matrix{<:Real})
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

    res = Vector{Matrix{Float64}}()
    
    for p in pint
        push!(res, reduce(hcat, p[1])')
    end

    return res    
end

function polyintersect(verts1::Vector{<:Vector{<:Real}}, verts2::Vector{<:Vector{<:Real}})
    return polyintersect(Matrix{Float64}(reduce(hcat, verts1)'), Matrix{Float64}(reduce(hcat, verts2)'))
end

end # module