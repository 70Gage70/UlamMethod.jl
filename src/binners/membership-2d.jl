function _membership2d(data::Matrix{<:Real}, bins::Bins{2, CRS}) where {CRS}
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

function _membership2d(traj::Trajectories{2}, bins::Bins{2, CRS}) where {CRS}    
    n_points = size(traj.x0, 2)
    membs = _membership2d(hcat(traj.x0, traj.xT), bins)
    return (membs[1:n_points], membs[n_points+1:end])
end
