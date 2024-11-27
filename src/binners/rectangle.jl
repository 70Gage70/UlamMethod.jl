"""
    struct RectangleBinner{M, CRS}

Bin a two dimensional `Polygon` with a tight covering of rectangles. The rectangles
are chosen to be as close as possible to squares and so the final number of bins
may be slightly different than the number requested.

### Fields

- `boundary`: A [`Boundary`](@ref) object.
- `bins`: A [`Bins`](@ref) object.
- `idx2pos`: A vector such that `idx2pos[i]` gives the position (in bins) of the bin initially \
(before removing dataless and disconnected bins) labelled `i`. `idx2pos[i] == nothing` if this bin was removed.

### Constructor 

    RectangleBinner(nbins, boundary; hardclip = true)

`nbins` can be:

- An `Integer`, in which case a `RectangleBinner` will be constructed attempting to distribute boxes proportionally \
across `x` and y` such that the total number of boxes is as close as possible to `nbins`.
- An `NTuple{2, Integer}`, in which case dimension `i` receives `nbins[i]` bins.
- An `NTuple{2, AbstractVector}`, in which case `nbins[i]` defines the edges of the bins in dimension `i`.
"""
struct RectangleBinner{M, CRS} <: BinningAlgorithm{2}
    boundary::Boundary{2, M, CRS}
    bins::Bins{2, M, CRS}
    idx2pos::Vector{Union{Int64, Nothing}}
end

function RectangleBinner(
    nbins::Union{Integer, NTuple{2, Integer}, NTuple{2, AbstractVector}}, 
    boundary::Boundary{2, M, CRS}; 
    hardclip::Bool = true) where {M, CRS}
    bbox = boundingbox(boundary.boundary)
    bbox_W = coords(bbox.max).x - coords(bbox.min).x
    bbox_L = coords(bbox.max).y - coords(bbox.min).y

    if nbins isa Integer
        nbins = hardclip ? ceil(Int64, nbins*area(bbox)/area(boundary.boundary)) : nbins
        n_x = sqrt(bbox_W*nbins/bbox_L) |> x -> round(Int64, x)
        n_y = sqrt(bbox_L*nbins/bbox_W) |> x -> round(Int64, x)
        grid = CartesianGrid(bbox.min, bbox.max, dims = (n_x, n_y))
    elseif nbins isa NTuple{2, Integer}
        grid = CartesianGrid(bbox.min, bbox.max, dims = nbins)
    elseif nbins isa NTuple{2, AbstractVector}
        @argcheck all(issorted.(nbins))
        grid = RectilinearGrid(nbins[1], nbins[2])
    end
    
    bins = Polytope{2, 𝔼{2}, CRS}[]

    for bin_ in elements(grid)
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

    return RectangleBinner(boundary, Bins(bins), Vector{Union{Int64, Nothing}}(1:length(bins)))
end 

membership(data::Matrix{<:Real}, binner::RectangleBinner{M, CRS}) where {M, CRS} = _membership2d(data, binner.bins)
membership(traj::Trajectories{2}, binner::RectangleBinner{M, CRS}) where {M, CRS} = _membership2d(traj, binner.bins)