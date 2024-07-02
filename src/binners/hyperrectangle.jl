"""
    struct HyperRectangleBinner{Dim, CRS}

Bin `Dim`-dimensional (where `Dim ≥ 3`) dataset based on a tight covering of identical HyperRectangles.

The hyperrectangles are chosen to be as close as possible to hypercubes.

The `HyperRectangleBinner` does not have hard clipping functionality. However, bins with no data are 
still removed in the main Ulam's method calculation.

The final number of bins may be different (larger) than the number requested.

### Fields 

- `boundary`: A [`Boundary`](@ref) object.
- `bins`: A [`Bins`](@ref) object.
- `ranges`: A vector such that `ranges[i]` is the `AbstractRange` giving the gridpoints of the bins along dimension `i`.

### Constructor

`HyperRectangleBinner(nbins, boundary)`
"""
struct HyperRectangleBinner{Dim, CRS, R<:AbstractRange} <: BinningAlgorithm{Dim}
    boundary::Boundary{Dim, CRS}
    bins::Bins{Dim, CRS}
    ranges::Vector{R}
end

function HyperRectangleBinner(nbins::Int64, boundary::Boundary{Dim, CRS}) where {Dim, CRS}
    @argcheck Dim >= 3 "In dimensions 1 and 2, LineBinner and RectangleBinner are preferred (respectively)."

    x_min = boundary.boundary.min
    x_max = boundary.boundary.max
    side_lengths = x_max .- x_min
    n_per_dim = [1 + ceil(Int64, side_length*(nbins/prod(side_lengths))^(1/Dim)) for side_length in side_lengths]
    # add 1 to account for the fact that n gridpoints makes n - 1 boxes

    # place vertices on grid 
    ranges = [range(x_min[i], x_max[i], length = n_per_dim[i]) for i = 1:Dim]

    verts = Iterators.product(ranges...) |> collect
    idx = Iterators.product(axes(verts)...) |> collect

    bins = Polytope{Dim, Dim, CRS}[]

    for i in idx
        if all((i .+ 1) .<= length.(ranges))
            j = i .+ 1
            bin_min, bin_max = verts[i...], verts[j...]
            hc = HyperRectangle(bin_min, bin_max)
            push!(bins, hc)
        end
    end

    return HyperRectangleBinner(boundary, Bins(bins), ranges)
end


function membership(data::Matrix{<:Real}, binner::HyperRectangleBinner{Dim, CRS, R}) where {Dim, CRS, R}
    @argcheck size(data, 1) == Dim

    ranges = binner.ranges
    verts = Iterators.product([r[1:end-1] for r in ranges]...) |> collect # N + 1 vertices is N bins so use [1:end-1]
    membs = Union{Int64, Nothing}[]

    for i = 1:size(data, 2)
        pt = data[:,i]
        dim_idx = zeros(Int64, Dim)
        in_nirv = false
        for j = 1:Dim
            if !(first(ranges[j]) ≤ pt[j] ≤ last(ranges[j]))
                in_nirv = true
                continue
            else
                dim_idx[j] = ceil(Int64, (pt[j] - first(ranges[j]))/step(ranges[j]))
            end
        end

        if in_nirv
            push!(membs, nothing)
        else
            push!(membs, LinearIndices(verts)[dim_idx...])
        end
    end

    return membs
end

function membership(traj::Trajectories{Dim}, binner::HyperRectangleBinner{Dim, CRS, R}) where {Dim, CRS, R}    
    n_points = size(traj.x0, 2)
    membs = membership(hcat(traj.x0, traj.xT), binner)
    return (membs[1:n_points], membs[n_points+1:end])
end