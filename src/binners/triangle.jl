"""
    struct TriangleBinner{CRS}

Bin a two dimensional `Polygon` with a covering of equilateral triangles. 

The final number of bins may be slightly different than the number requested.

### Fields

- `boundary`: A [`Boundary`](@ref) object.
- `bins`: A [`Bins`](@ref) object.
- `idx2pos`: A vector such that `idx2pos[i]` gives the position (in bins) of the bin initially \
(before removing dataless and disconnected bins) labelled `i`. `idx2pos[i] == nothing` if this bin was removed.

### Constructor 

`TriangleBinner(nbins, boundary; hardclip = true)`
"""
struct TriangleBinner{CRS} <: BinningAlgorithm{2}
    boundary::Boundary{2, CRS}
    bins::Bins{2, CRS}
    idx2pos::Vector{Union{Int64, Nothing}}
end

function TriangleBinner(nbins::Int64, boundary::Boundary{2, CRS}; hardclip::Bool = true) where {CRS}
    bbox = Meshes.boundingbox(boundary.boundary)
    xmin, xmax, ymin, ymax = coords(bbox.min).x.val, coords(bbox.max).x.val, coords(bbox.min).y.val, coords(bbox.max).y.val
    W = xmax - xmin
    L = ymax - ymin

    nbins = hardclip ? ceil(Int64, (1/6)*nbins*area(bbox)/area(boundary.boundary)) : nbins
    
    α = L/W
    β = sqrt(1 + 24*sqrt(3)*α*nbins)
    n_x, n_y = ceil(Int64, (β - 1)/(4*sqrt(3)*α)), ceil(Int64, (β + 1)/6)
    R = max(W/(sqrt(3)*n_x), 2*L/(3*n_y - 1))
    r = sqrt(3)*R/2
    
    cs = []
    
    for y = 1:n_y
        x_left = isodd(y) ? 0.0 : -r
        n_x_corrected = isodd(y) ? n_x : n_x + 1
        for x = 1:n_x_corrected
            push!(cs, [x_left + (x - 1)*2*r, (y - 1)*3*R/2])
        end
    end 
    
    triangles = []
    
    for c in cs, t = 1:6
        push!(triangles, [
            (c[1], c[2]), 
            (c[1] + R*sin(π*t/3), c[2] + R*cos(π*t/3)), 
            (c[1] + R*sin(π*(t + 1)/3), c[2] + R*cos(π*(t + 1)/3))
            ]
        )
    end
    
    triangles = [Ngon(hex...) for hex in triangles]
    
    c_box = ((xmin + xmax)/2, (ymin + ymax)/2)
    c_hex = centroid(boundingbox(triangles)) |> p -> (coords(p).x.val, coords(p).y.val)
    triangles = Translate(c_box .- c_hex)(triangles)

    # points near but not exactly equal to zero cause problems for Meshes.intersect, so round them to zero
    near_zero_p(p) = ([coords(p).x.val, coords(p).y.val] .|> x -> (y -> abs(y) < 10^(-15) ? 0.0 : y).(x)) |> z -> Meshes.Point(z...)
    triangles = [Triangle(near_zero_p.(vertices(tri))...) for tri in triangles]

    bins = Polytope{2, 2, CRS}[]
    for bin_ in triangles
        isect = hardclip ? intersect(bin_, boundary.boundary) : intersects(bin_, boundary.boundary) ? bin_ : nothing
        if isect isa PolyArea
            # if the intersection is a PolyArea (i.e. a polygon possibly with holes), we convert it to an Ngon with the outer ring
            verts = isect.rings[1].vertices
            if length(verts) >= 3
                push!(bins, Ngon(verts...))
            end
        elseif isect isa Ngon
            push!(bins, isect)
        end
    end

    return TriangleBinner(boundary, Bins(bins), Vector{Union{Int64, Nothing}}(1:length(bins)))
end 

membership(data::Matrix{<:Real}, binner::TriangleBinner{CRS}) where {CRS} = _membership2d(data, binner.bins)
membership(traj::Trajectories{2}, binner::TriangleBinner{CRS}) where {CRS} = _membership2d(traj, binner.bins)