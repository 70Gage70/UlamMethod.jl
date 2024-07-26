"""
    struct LineBinner{M, CRS}

Bin a one dimensional `Segment` (line segement) with `nbins` equally-spaced bins.

### Fields

- `boundary`: A [`Boundary`](@ref) object.
- `bins`: A [`Bins`](@ref) object.
- `idx2pos`: A vector such that `idx2pos[i]` gives the position (in bins) of the bin initially \
(before removing dataless and disconnected bins) labelled `i`. `idx2pos[i] == nothing` if this bin was removed.

### Constructor 

    LineBinner(nbins, boundary; hardclip = true)
"""
struct LineBinner{M, CRS} <: BinningAlgorithm{1}
    boundary::Boundary{1, M, CRS}
    bins::Bins{1, M, CRS}
    idx2pos::Vector{Union{Int64, Nothing}}
end

function LineBinner(nbins::Int64, boundary::Boundary{1, M, CRS}; hardclip::Bool = true) where {M, CRS}
    bbox = Meshes.boundingbox(boundary.boundary)
    grid = CartesianGrid(bbox.min, bbox.max, dims = (nbins, ))

    bins = Polytope{1, 𝔼{1}, CRS}[]
    for bin_ in elements(grid)
        isect = hardclip ? intersect(bin_, boundary.boundary) : intersects(bin_, boundary.boundary) ? bin_ : nothing
        if isect isa Segment
            push!(bins, isect)
        end
    end

    return LineBinner(boundary, Bins(bins), Vector{Union{Int64, Nothing}}(1:length(bins)))
end 

membership(data::Matrix{<:Real}, binner::LineBinner{M, CRS}) where {M, CRS} = _membership1d(data, binner.bins)
membership(traj::Trajectories{1}, binner::LineBinner{M, CRS}) where {M, CRS} = _membership1d(traj, binner.bins)