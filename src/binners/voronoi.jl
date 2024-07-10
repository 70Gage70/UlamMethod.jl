"""
    struct VoronoiBinner{Dim, CRS}

Bin `Dim`-dimensional (where `Dim ∈ (1, 2)`) dataset based on a Voronoi tesellation of points \
generated by a k-means clustering of initial trajectory data (i.e. `Trajectories.x0`).

### Fields 

- `boundary`: A [`Boundary`](@ref) object.
- `bins`: A [`Bins`](@ref) object.
- `idx2pos`: A vector such that `idx2pos[i]` gives the position (in bins) of the bin initially \
(before removing dataless and disconnected bins) labelled `i`. `idx2pos[i] == nothing` if this bin was removed.

### Constructor

`VoronoiBinner(nbins, boundary, traj; hardclip = true)`
"""
struct VoronoiBinner{Dim, CRS} <: BinningAlgorithm{Dim}
    boundary::Boundary{Dim, CRS}
    bins::Bins{Dim, CRS}
    idx2pos::Vector{Union{Int64, Nothing}}
end

function VoronoiBinner(
    nbins::Int64, 
    boundary::Boundary{Dim, CRS}, 
    traj::Trajectories{Dim}; 
    hardclip::Bool = true) where {Dim, CRS}

    @argcheck Dim in [1, 2]

    nbins = nbins
    pts = [points(boundary) ;; traj.x0]

    centers = kmeans(Yinyang(), pts, nbins).centers |> x -> [Point(x[:, i]...) for i = 1:size(x, 2)]

    if Dim == 1
        grid = RectilinearGrid([coords(c).x for c in centers])
    elseif Dim == 2
        grid = tesselate(centers, VoronoiTesselation())
    end

    bins = Polytope{Dim, Dim, CRS}[]

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

    return VoronoiBinner(boundary, Bins(bins), Vector{Union{Int64, Nothing}}(1:length(bins)))
end 

membership(data::Matrix{<:Real}, binner::VoronoiBinner{1, CRS}) where {CRS} = _membership1d(data, binner.bins)
membership(traj::Trajectories{1}, binner::VoronoiBinner{1, CRS}) where {CRS} = _membership1d(traj, binner.bins)

membership(data::Matrix{<:Real}, binner::VoronoiBinner{2, CRS}) where {CRS} = _membership2d(data, binner.bins)
membership(traj::Trajectories{2}, binner::VoronoiBinner{2, CRS}) where {CRS} = _membership2d(traj, binner.bins)