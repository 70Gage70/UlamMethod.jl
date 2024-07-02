"""
    struct LineBinner{CRS}

Bin a one dimensional `Segment` (line segement) with `nbins` equally-spaced bins.

### Fields

- `boundary`: A [`Boundary`](@ref) object.
- `bins`: A [`Bins`](@ref) object.

### Constructor 

`LineBinner(nbins, boundary; hardclip = true)`
"""
struct LineBinner{CRS} <: BinningAlgorithm{1}
    boundary::Boundary{1, CRS}
    bins::Bins{1, CRS}
end

function LineBinner(nbins::Int64, boundary::Boundary{1, CRS}; hardclip::Bool = true) where {CRS}
    bbox = Meshes.boundingbox(boundary.boundary)
    grid = CartesianGrid(bbox.min, bbox.max, dims = (nbins, ))

    bins = Polytope{1, 1, CRS}[]
    for bin_ in elements(grid)
        isect = hardclip ? intersect(bin_, boundary.boundary) : intersects(bin_, boundary.boundary) ? bin_ : nothing
        if isect isa Segment
            push!(bins, isect)
        end
    end

    return LineBinner(boundary, Bins(bins))
end 

membership(data::Matrix{<:Real}, binner::LineBinner{CRS}) where {CRS} = _membership1d(data, binner.bins)
membership(traj::Trajectories{1}, binner::LineBinner{CRS}) where {CRS} = _membership1d(traj, binner.bins)