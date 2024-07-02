"""
    struct HexagonBinner{CRS}

Bin a two dimensional `Polygon` with a covering of regular hexagons. 

The final number of bins may be slightly different than the number requested.

### Fields

- `boundary`: A [`Boundary`](@ref) object.
- `bins`: A [`Bins`](@ref) object.

### Constructor 

`HexagonBinner(nbins, boundary; hardclip = true)`
"""
struct HexagonBinner{CRS} <: BinningAlgorithm{2}
    boundary::Boundary{2, CRS}
    bins::Bins{2, CRS}
end

function HexagonBinner(nbins::Int64, boundary::Boundary{2, CRS}; hardclip::Bool = true) where {CRS}
    bbox = Meshes.boundingbox(boundary.boundary)
    xmin, xmax, ymin, ymax = coords(bbox.min).x.val, coords(bbox.max).x.val, coords(bbox.min).y.val, coords(bbox.max).y.val
    W = xmax - xmin
    L = ymax - ymin

    nbins = hardclip ? ceil(Int64, nbins*area(bbox)/area(boundary.boundary)) : nbins
    
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
    
    hexagons = []
    
    for c in cs
        push!(hexagons, [(c[1] + R*sin(π*k/3), c[2] + R*cos(π*k/3)) for k = 0:5])
    end
    
    hexagons = [Ngon(hex...) for hex in hexagons]
    
    c_box = ((xmin + xmax)/2, (ymin + ymax)/2)
    c_hex = centroid(boundingbox(hexagons)) |> p -> (coords(p).x.val, coords(p).y.val)
    hexagons = Translate(c_box .- c_hex)(hexagons)

    bins = Polytope{2, 2, CRS}[]
    for bin_ in hexagons
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

    return HexagonBinner(boundary, Bins(bins))
end 

membership(data::Matrix{<:Real}, binner::HexagonBinner{CRS}) where {CRS} = _membership2d(data, binner.bins)
membership(traj::Trajectories{2}, binner::HexagonBinner{CRS}) where {CRS} = _membership2d(traj, binner.bins)