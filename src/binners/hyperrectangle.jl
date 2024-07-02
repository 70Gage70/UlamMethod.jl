"""
    struct HyperRectangleBinner{Dim}

Bin `Dim`-dimensional (where `Dim ≥ 3`) dataset based on a tight covering of identical HyperRectangles.

The hyperrectangles are chosen to be as close as possible to hypercubes.

The `HyperRectangleBinner` does not have hard clipping functionality. However, bins with no data are 
still removed in the main Ulam's method calculation.

The final number of bins may be different (larger) than the number requested.

### Fields 

- `nbins`: The number of bins requested.
- `boundary`: The `Boundary` defining the domain.
- `ranges`: A vector such that `ranges[i]` is the `AbstractRange` giving the gridpoints of the bins along dimension `i`.

### Constructor

`HyperRectangleBinner(nbins, boundary)`
"""
struct HyperRectangleBinner{Dim, CRS, R<:AbstractRange} <: BinningAlgorithm{Dim}
    nbins::Int64
    boundary::Boundary{Dim, CRS}
    ranges::Vector{R}
    
    function HyperRectangleBinner(nbins::Int64, boundary::Boundary{Dim, CRS}) where {Dim, CRS}
        @argcheck Dim >= 3 "In dimensions 1 and 2, LineBinner and RectangleBinner are preferred (respectively)."

        x_min = boundary.boundary.min
        x_max = boundary.boundary.max
        side_lengths = x_max .- x_min
        n_per_dim = [1 + ceil(Int64, side_length*(nbins/prod(side_lengths))^(1/Dim)) for side_length in side_lengths]
        # add 1 to account for the fact that n gridpoints makes n - 1 boxes
    
        # place vertices on grid 
        ranges = [range(x_min[i], x_max[i], length = n_per_dim[i]) for i = 1:Dim]

        new{Dim, CRS, eltype(ranges)}(nbins, boundary, ranges)
    end
end

function _bin(boundary::Boundary{Dim, CRS}, binner::HyperRectangleBinner{Dim}) where {Dim, CRS}
    verts = Iterators.product(binner.ranges...) |> collect
    idx = Iterators.product(axes(verts)...) |> collect

    bins = Polytope{Dim, CRS}[]

    for i in idx
        if all((i .+ 1) .<= length.(binner.ranges))
            j = i .+ 1
            bin_min, bin_max = verts[i...], verts[j...]
            hc = HyperRectangle(bin_min, bin_max)
            push!(bins, hc)
        end
    end

    return Bins(bins)
end

### ≥3D Algorithm

# function _membershipND(data::Matrix{<:Real}, bins::Bins{K, Dim, CRS}) where {K, Dim, CRS}
#     @argcheck size(data, 1) == Dim

#     membs = Int64[]

#     # for bin in bins.bins


#     return membs
# end