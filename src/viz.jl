@doc """
    viz(object; kwargs...)

This function wraps `Meshes.viz` for `UlamMethod` objects for plotting.

### Methods

    viz(boundary; kwargs...)

Plot a boundary.

    viz(bins; kwargs...)

Plot the bins.

    viz(binner; kwargs...)

Plot `binner.bins`.

    viz(data; kwargs...)

Plot a `2 x N` matrix of points.

    viz(traj; kwargs...)

Plot `traj.x0` and `traj.xT` on the same graph.

    viz(ulam; bins_labels = true, bins_labels_size = 16)

Plot the bins, boundary and a heatmap of the stationary distribution of `P_closed(ulam)`.

Optionally, superimpose bin labels with `bin_labels` and adjust their size with `bins_labels_size`.

    viz(traj, ulam; bins_labels = true, bins_labels_size = 16)

Plot a heatmap of the counts in each bin of both `traj.x0` and `traj.xT`.

Optionally, superimpose bin labels with `bin_labels` and adjust their size with `bins_labels_size`.

---

### The docstring for `Meshes.viz` is provided below.

---

$(@doc(Meshes.viz))
"""
function viz end

@doc """
    viz!(object; kwargs...)

This function wraps `Meshes.viz!` for `UlamMethod` objects for plotting.

### Methods

    viz!(boundary; kwargs...)

Plot a boundary.

    viz!(bins; kwargs...)

Plot the bins.

    viz!(binner; kwargs...)

Plot `binner.bins`.

    viz!(data; kwargs...)

Plot a `2 x N` matrix of points.

---

### The docstring for `Meshes.viz!` is provided below.

---

$(@doc(Meshes.viz!))
"""
function viz! end