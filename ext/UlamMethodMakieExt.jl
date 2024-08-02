module UlamMethodMakieExt

using UlamMethod
using ArgCheck
using Makie
using LinearAlgebra: eigvecs, normalize 
import Meshes

### BOUNDARY
UlamMethod.viz(boundary::Boundary; kwargs...) = Meshes.viz(boundary.boundary; merge((showsegments = true,), kwargs)...)
UlamMethod.viz!(boundary::Boundary; kwargs...) = Meshes.viz!(boundary.boundary; merge((showsegments = true,), kwargs)...)

### BINS
UlamMethod.viz(bins::Bins; kwargs...) = Meshes.viz(bins.bins; merge((showsegments = true,), kwargs)...)
UlamMethod.viz!(bins::Bins; kwargs...) = Meshes.viz!(bins.bins; merge((showsegments = true,), kwargs)...)

### BINNER
UlamMethod.viz(binner::BinningAlgorithm; kwargs...) = UlamMethod.viz(binner.bins; kwargs...)
UlamMethod.viz!(binner::BinningAlgorithm; kwargs...) = UlamMethod.viz!(binner.bins; kwargs...)

### RESULT
function UlamMethod.viz(
    ulam::UlamResult{Dim, M, CRS}; 
    bins_labels::Bool = true,
    bins_labels_size::Integer = 16) where {Dim, M, CRS}

    @argcheck Dim == 2 "Plotting is only available in 2D"
    @argcheck labels_size > 0

    fig = Figure()
    ax = Axis(fig[1, 1], 
        xgridvisible = false,
        ygridvisible = false,
        xlabel = "x",
        ylabel = "y")

    ### BINS
    pi_stat = abs.(eigvecs(P_closed(ulam)')[:,end])[1:end-1] |> x -> normalize(x, 1)
    UlamMethod.viz!(ulam.binner, color = pi_stat, segmentcolor = :white)

    ### BINS_DIS
    length(ulam.bins_dis.bins) > 0 && UlamMethod.viz!(ulam.bins_dis, color = :black, segmentcolor = :white)

    ### BOUNDARY
    # pb = points(ulam.binner.boundary)
    # ls = [(pb[1, i], pb[2, i]) for i = 1:size(pb, 2)] |> x -> [(x[i], x[i + 1]) for i = 1:length(x)-1]
    # ls = [ls ; ((pb[1, end], pb[2, end]), (pb[1, end], pb[2, end]))]
    UlamMethod.viz!(ulam.binner.boundary, segmentsize = 5, segmentcolor = :black, color = :white, alpha = 0.0)
    
    ### COLORBAR
    Colorbar(fig[1, 2], colorrange = (0, maximum(pi_stat)), label = "Stationary distribution weight")

    if bins_labels
        centers = [(c.coords.x.val, c.coords.y.val) for c in Meshes.centroid.(bins(ulam).bins)]
        text!(centers, text = string.(1:length(centers)), 
            align = (:center, :center), 
            color = :white,
            fontsize = bins_labels_size)
    end

    fig
end


end # module