"""
    struct HyperRectangle{Dim, M, CRS}

An extention of `Meshes` to dimensions greater than 3.

The extension is minimal and implements no methods and does not interface with CRS.

### Fields 

- `min`: The coordinates of the minimum of the HyperRectangle as an `NTuple`.
- `max`: The coordinates of the maximum of the HyperRectangle as an `NTuple`.
- `vertices`: `[min, max]`, required by `Meshes`.

### Constructor

`HyperRectangle(min, max)`
"""
struct HyperRectangle{Dim, M, CRS} <: Polytope{Dim, M, CRS}
    min::NTuple{Dim, Float64}
    max::NTuple{Dim, Float64}
    vertices::Vector{NTuple{Dim, Float64}}

    function HyperRectangle(min::NTuple{Dim, Real}, max::NTuple{Dim, Real}) where {Dim}
        b = Box(tuple(zeros(Dim)...), tuple(zeros(Dim)...)) |> typeof
        return new{Dim, b.parameters[1], b.parameters[2]}(min, max, [min, max])
    end
end

"""
    struct HyperRectangleBinner{Dim, M, CRS}

Bin `Dim`-dimensional (where `Dim ≥ 3`) dataset based on a tight covering of identical HyperRectangles.

The hyperrectangles are chosen to be as close as possible to hypercubes.

The `HyperRectangleBinner` does not have hard clipping functionality. However, bins with no data are 
still removed in the main Ulam's method calculation.

The final number of bins may be different (larger) than the number requested.

### Fields 

- `boundary`: A [`Boundary`](@ref) object.
- `bins`: A [`Bins`](@ref) object.
- `idx2pos`: A vector such that `idx2pos[i]` gives the position (in bins) of the bin initially \
(before removing dataless and disconnected bins) labelled `i`. `idx2pos[i] == nothing` if this bin was removed.
- `ranges`: A vector such that `ranges[i]` is the `AbstractRange` giving the gridpoints of the bins along dimension `i`.

### Constructor

`HyperRectangleBinner(nbins, boundary)`
"""
struct HyperRectangleBinner{Dim, M, CRS, R<:AbstractRange} <: BinningAlgorithm{Dim}
    boundary::Boundary{Dim, M, CRS}
    bins::Bins{Dim, M, CRS}
    idx2pos::Vector{Union{Int64, Nothing}}
    ranges::Vector{R}
end

function HyperRectangleBinner(nbins::Int64, boundary::Boundary{Dim, M, CRS}) where {Dim, M, CRS}
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

    bins = Polytope{Dim, 𝔼{Dim}, CRS}[]

    for i in idx
        if all((i .+ 1) .<= length.(ranges))
            j = i .+ 1
            bin_min, bin_max = verts[i...], verts[j...]
            hc = HyperRectangle(bin_min, bin_max)
            push!(bins, hc)
        end
    end

    return HyperRectangleBinner(boundary, Bins(bins), Vector{Union{Int64, Nothing}}(1:length(bins)), ranges)
end


function membership(data::Matrix{<:Real}, binner::HyperRectangleBinner{Dim, M, CRS, R}) where {Dim, M, CRS, R}
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

    # this entire calculation assumes that the bins are arranged as they are when `HyperRectangleBinner` is 
    # initially constructed; when bins are removed, the structure changes, so we have to find where the bins 
    # actually are using `idx2pos`
    return [isnothing(m) ? nothing : binner.idx2pos[m] for m in membs]
end

function membership(traj::Trajectories{Dim}, binner::HyperRectangleBinner{Dim, M, CRS, R}) where {Dim, M, CRS, R}    
    return (membership(traj.x0, binner), membership(traj.xT, binner))
end