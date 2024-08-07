module UlamMethodMakieExt

using UlamMethod
using ArgCheck
using Makie
using LinearAlgebra: eigvecs, normalize 
import Meshes
using PrecompileTools: @compile_workload 

### BOUNDARY
UlamMethod.viz(boundary::Boundary; kwargs...) = Meshes.viz(boundary.boundary; 
    merge((showsegments = true, segmentsize = 5, segmentcolor = :black, color = :white, alpha = 0.0), kwargs)...)
UlamMethod.viz!(boundary::Boundary; kwargs...) = Meshes.viz!(boundary.boundary; 
    merge((showsegments = true, segmentsize = 5, segmentcolor = :black, color = :white, alpha = 0.0), kwargs)...)

### BINS
UlamMethod.viz(bins::Bins; kwargs...) = Meshes.viz(bins.bins; merge((showsegments = true,), kwargs)...)
UlamMethod.viz!(bins::Bins; kwargs...) = Meshes.viz!(bins.bins; merge((showsegments = true,), kwargs)...)

### BINNER
UlamMethod.viz(binner::BinningAlgorithm; kwargs...) = UlamMethod.viz(binner.bins; kwargs...)
UlamMethod.viz!(binner::BinningAlgorithm; kwargs...) = UlamMethod.viz!(binner.bins; kwargs...)

### DATA
function UlamMethod.viz(data::Matrix{<:Real}; kwargs...)
    @argcheck size(data, 1) == 2 "Plotting is only available in 2 dimensions."
    ps = Meshes.PointSet(Meshes.Point(data[:,i]...) for i = 1:size(data, 2))

    Meshes.viz(ps; kwargs...)
end

function UlamMethod.viz!(data::Matrix{<:Real}; kwargs...)
    @argcheck size(data, 1) == 2 "Plotting is only available in 2 dimensions."
    ps = Meshes.PointSet(Meshes.Point(data[:,i]...) for i = 1:size(data, 2))

    Meshes.viz!(ps; kwargs...)
end

### TRAJECTORY
function UlamMethod.viz(traj::Trajectories{Dim}) where {Dim}
    @argcheck Dim == 2 "Plotting is only available in 2D"

    fig = Figure()
    ax = Axis(fig[1, 1], 
        xlabel = "x",
        ylabel = "y")
    x0 = traj.x0 |> x -> Meshes.PointSet(Meshes.Point(x[:,i]...) for i = 1:size(x, 2))
    xT = traj.xT |> x -> Meshes.PointSet(Meshes.Point(x[:,i]...) for i = 1:size(x, 2))

    x0_plot = Meshes.viz!(x0, color = :blue)
    xT_plot = Meshes.viz!(xT, color = :red)

    Legend(fig[1, 2], [x0_plot, xT_plot], ["x0", "xT"])

    fig
end

### RESULT
function UlamMethod.viz(
    ulam::UlamResult{Dim, M, CRS}; 
    bins_labels::Bool = true,
    bins_labels_size::Integer = 16) where {Dim, M, CRS}

    @argcheck Dim == 2 "Plotting is only available in 2D"
    @argcheck bins_labels_size > 0

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
    UlamMethod.viz!(ulam.binner.boundary)
    
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

function UlamMethod.viz(
    traj::Trajectories{Dim},
    ulam::UlamResult{Dim, M, CRS}; 
    bins_labels::Bool = true,
    bins_labels_size::Integer = 16) where {Dim, M, CRS}

    @argcheck Dim == 2 "Plotting is only available in 2D"
    @argcheck bins_labels_size > 0

    fig = Figure(size = (1000, 500))

    ### X0
    ax1 = Axis(fig[1, 1],
        aspect = DataAspect(), 
        xgridvisible = false,
        ygridvisible = false,
        xlabel = "x",
        ylabel = "y")

    # BINS
    cts = counts(traj.x0, ulam)[1:end-1]
    UlamMethod.viz!(ulam.binner, color = cts, segmentcolor = :white)

    # BINS_DIS
    length(ulam.bins_dis.bins) > 0 && UlamMethod.viz!(ulam.bins_dis, color = :black, segmentcolor = :white)

    # BOUNDARY
    UlamMethod.viz!(ulam.binner.boundary)
    
    # COLORBAR
    Colorbar(fig[2, 1], colorrange = (0, maximum(cts)), label = "x0 Counts", vertical = false)

    ### XT
    ax2 = Axis(fig[1, 2],
        aspect = DataAspect(), 
        xgridvisible = false,
        ygridvisible = false,
        xlabel = "x",
        ylabel = "y")

    # BINS
    cts = counts(traj.xT, ulam)[1:end-1]
    UlamMethod.viz!(ulam.binner, color = cts, segmentcolor = :white)

    # BINS_DIS
    length(ulam.bins_dis.bins) > 0 && UlamMethod.viz!(ulam.bins_dis, color = :black, segmentcolor = :white)

    # BOUNDARY
    UlamMethod.viz!(ulam.binner.boundary)
    
    # COLORBAR
    Colorbar(fig[2, 2], colorrange = (0, maximum(cts)), label = "xT Counts", vertical = false)


    if bins_labels
        centers = [(c.coords.x.val, c.coords.y.val) for c in Meshes.centroid.(bins(ulam).bins)]
        text!(ax1, centers, text = string.(1:length(centers)), 
            align = (:center, :center), 
            color = :white,
            fontsize = bins_labels_size)
        text!(ax2, centers, text = string.(1:length(centers)), 
        align = (:center, :center), 
        color = :white,
        fontsize = bins_labels_size)
    end

    fig
end

@compile_workload begin
    import Random
    Random.seed!(1234)

    traj = Trajectories(2, 1000)
    ulam = ulam_method(traj, 200)

    viz(traj)
    viz(ulam)
end

end # module