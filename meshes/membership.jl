### 1D Algorithm
function _membership(data::Matrix{<:Real}, bins::Bins{K, 1, CRS}) where {K, CRS}   
    @argcheck size(data, 1) == 1 
    memb = Union{Nothing, Int64}[]

    for x0 in data
        memb_i = nothing
        for i = 1:length(bins.bins)
            bin_ = bins.bins[i].vertices .|> x -> coords(x).x.val
            if bin_[1] < x0 < bin_[2]
                memb_i = i
                continue
            end
        end
        push!(memb, memb_i)
    end

    return memb
end

function _membership(traj::Trajectories{1}, bins::Bins{K, 1, CRS}) where {K, CRS}    
    return (_membership(traj.x0, bins), _membership(traj.xT, bins))
end

### 2D Algorithm
function _membership(data::Matrix{<:Real}, bins::Bins{K, 2, CRS}) where {K, CRS}
    @argcheck size(data, 1) == 2

    nodes = []
    edges = []
    
    for i = 1:length(bins.bins)
        bin_i = bins.bins[i]
        push!(nodes, collect(coords.(bin_i.vertices) .|> x -> [x.x.val, x.y.val]))
        node_idx = i == 1 ? 0 : edges[end][end][1]
        n_nodes = length(nodes[i])
        push!(edges, [[node_idx + k, k == n_nodes ? node_idx + 1 : node_idx + k + 1, i] for k = 1:n_nodes])
    end
    
    nodes = vcat(stack.(nodes, dims = 1)...)
    edges = vcat(stack.(edges, dims = 1)...)
    
    n_points = size(data, 2)
    ip2 = inpoly2(permutedims(data), nodes, edges) 
    # returns stats[:, 1:2, area], where : has the index of the points and 1:2 are [inside, onboundary]
    # therefore, stats[:, 1, k] is a bitvector of point membership to the kth polygon
    # therefore, stats[n, 1, :] is bitvector of polygon membership to the nth point
    
    membs = [
        begin 
            # the polygons don't overlap, so each point only belongs to one polygon and we can use `findfirst` to find that polygon's index
            ff = findfirst(ip2[pt,1,:]); 
            isnothing(ff) ? nothing : ip2[pt,2,ff] ? nothing : ff # if a point is on an edge (ip2[pt, 2, ff] == true), we ignore it
        end 
        for pt = 1:n_points
    ]

    return membs
end

function _membership(traj::Trajectories{2}, bins::Bins{K, 2, CRS}) where {K, CRS}    
    n_points = size(traj.x0, 2)
    membs = _membership(hcat(traj.x0, traj.xT), bins)
    return (membs[1:n_points], membs[n_points+1:end])
end

"""
    membership(traj, bins)

Compute `(memb_x0, memb_xT)` where `memb_x0[i] = j` if `traj.x0[i]` is inside `bins.bins[j]` and `memb_x0[i] = nothing` if it is inside no bin. Similarly for `xT`.
"""
function membership(traj::Trajectories{Dim}, bins::Bins{K, Dim, CRS}) where {K, Dim, CRS}
    return _membership(traj, bins)
end