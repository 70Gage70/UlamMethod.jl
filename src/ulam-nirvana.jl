"""
A generic Ulam method.
Needs: x0, y0, xT, yT, polys
Needs: polys: a vector whose i'th entry is a vector of [x, y] coordinates of that polygon WITH CLOSING NODE
Needs: polys_centers: a vector whose i'th entry is the [x, y] coordinates of the center of the ith polygon. 

Written by Gage Bonner November 2022
"""
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

using PolygonInbounds:inpoly2
using Graphs
using LinearAlgebra

include("helpers.jl")

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

"""
Generates a list of nodes and edges suitable for use in inpoly2(). Here, polys is from the binners.
"""

function inpoly_preprocess(polys)
    n_poly = length(polys)
    nodes = zeros(sum(length(polys[i]) - 1 for i = 1:n_poly), 2) # x | y
    edges = zeros(Int64, sum(length(polys[i]) - 1 for i = 1:n_poly), 3) # x | y | cell index
    n_edge = 1
    edge_max = 0

    # have to construct nodes (defined by vertices of Voronoi polygons in polys) and edges (just integers linking vertices) for inpoly2
    for i = 1:n_poly
        for j = 1:length(polys[i])-1
            if j < length(polys[i])-1
                edges[n_edge,:] = [edge_max + j, edge_max + j + 1, i]
            else # have to close polygon
                edges[n_edge,:] = [edge_max + j, edge_max + 1, i]
                edge_max = edge_max + j
            end

            nodes[n_edge,:] = [polys[i][j][1], polys[i][j][2]]
            n_edge = n_edge + 1
        end
    end

    return Dict("nodes" => nodes, "edges" => edges)
end

"""
Get the largest strongly connected component. Construct the adjacency matrix corresponding to P_closed, use that to make the associated
directed graph, and use an scc algorithm from Graphs.jl
"""

function get_scc_inds(P_closed)
    Psize = length(P_closed[1,:])
    P_open_size = Psize - 1 # exclude nirvana when calculating scc, assumes one nirvana state in last entry
    
    Padj = Matrix{Int64}(undef, P_open_size, P_open_size) # the adjacency matrix defined by P
    for i = 1:P_open_size
        for j = 1:P_open_size
            if P_closed[i, j] == 0.0
                Padj[i,j] = 0 
            else
                Padj[i,j] = 1
            end
        end
    end
    
    Pgraph = SimpleDiGraph(Padj)
    scc = sort(strongly_connected_components(Pgraph), by = length) 
    largest_component = sort(scc[end])
    push!(largest_component, Psize) # put nirvana back at the end, should stay sorted

    return largest_component
end


"""
Processes the output of ulam_nirvana ("polys") for use with inpoly2. 
"""

function inpoly_preprocess_AB(polys)
    nodes = []
    edges = []

    poly_index = 1
    edge_index = 1
    edges_this = 0
    i = 1

    while i <= size(polys)[1] - 1
        line = polys[i,:]
        if polys[i + 2] == poly_index
            push!(nodes, line[2:3])
            push!(edges, [edge_index, edge_index + 1, poly_index])
            edges_this = edges_this + 1
            i = i + 1
        elseif (polys[i + 2] != poly_index) || (i == size(polys)[1] - 1) 
            # this is the connecting edge, note we did i + 2 on the last line since polys has the connecting last edge explicitly
            push!(nodes, line[2:3])
            push!(edges, [edge_index, edge_index - edges_this, poly_index])
            i = i + 2
            poly_index = poly_index + 1
            edges_this = 0
        end

        edge_index  = edge_index + 1
    end

    nodes = vecvec_to_mat(nodes)
    edges = vecvec_to_mat(edges)
    npolys = Integer(polys[end, 1])

    return Dict("nodes" => nodes, "edges" => edges, "npolys" => npolys)
end


"""
Finds the indices of ulam_polys which contain the points in centers. For finding the indices of A or B after ulam is already done.
"""

function getABinds(ulam_polys, centers)
    res = inpoly_preprocess_AB(ulam_polys)
    nodes_inpoly = res["nodes"]
    edges_inpoly = res["edges"]
    npolys = res["npolys"]

    inds = zeros(Int64, npolys)

    res_poly = inpoly2(centers, nodes_inpoly, edges_inpoly)

    for i = 1:npolys
        if (1 in res_poly[:,1,i]) # the i'th cell contains at least one point from A, therefore this cell belongs to A
            inds[i] = 1
        end
    end

    inds = findall(!iszero, inds)
    
    return inds
end


"""
The Ulam method. 

polys should be a vector of vectors; that is, poly[i] is a vector whose j'th entry is the jth vertex of that polygon. Must be closed!
polys_centers is also a vector of vectors; polys_centers[i] is equal to [x_center, y_center] 
"""

# ulam_nirvana(x0, y0, xT, yT, polys, polys_centers; sto_type = "data", source_centers = polys_centers[1])

function ulam_nirvana(traj::UlamTrajectories, domain::UlamDomain, polys::Vector{UlamPolygon})::UlamResult
    ##########################################################################################
    """
    Assign indices to obs/traj data based on given polys.
    Simultaneously purge polygons with no observational data in them.
    """
    ##########################################################################################
    # println("Assigning indices.")

    inds0 = inpoly([traj.x0 ;; traj.y0], PolyTable(polys))
    indsT = inpoly([traj.xT ;; traj.yT], PolyTable(polys))
    contains_data = []
    new_polys_count = 1


    polys_clean = polys[contains_data] # only keep the polygons that actually have data in them
    polys_centers_clean = polys_centers[contains_data]
    npolys = length(contains_data)

    if length(contains_data) == length(polys)
        clean_info = "All polygons contained data. "
    else
        clean_info = "Some polygons didn't contain data. A total of " * string(length(polys) - length(contains_data)) * " polygons were removed. "
    end

    ##########################################################################################
    """
    Construct P_closed.
    """
    ##########################################################################################

    Psize = npolys + 1 # extra row/col for nirvana
    P_closed = zeros(Psize, Psize)
    reinjection_counts = zeros(Int64, Psize)

    for i = 1:length(inds0)
        ind0 = inds0[i]
        indT = indsT[i]

        if ind0 != 0 && indT != 0 # transitions from interior to interior
            P_closed[ind0, indT] = P_closed[ind0, indT] + 1
        elseif ind0 != 0 && indT == 0 # transitions from interior to nirvana
            P_closed[ind0, Psize] = P_closed[ind0, Psize] + 1
        elseif ind0 == 0 && indT != 0 # transitions from nirvana to interior
            reinjection_counts[indT] = reinjection_counts[indT] + 1 # re-inject these later according to the desired type
        # else: transitions from nirvana to self; ignored            
        end
    end

    ####################################################################################################################################
    """
    Find the strongest connected component of P_closed.
    """
    ##########################################################################################
    # println("Extracting strongest connected component.")

    scc_inds = get_scc_inds(P_closed) # indices such that P_closed[scc_inds, scc_inds] is connected

    initial_size = Psize

    Psize = length(scc_inds)
    P_closed = P_closed[scc_inds, scc_inds]
    polys_clean_scc = splice!(polys_clean, scc_inds[1:end-1]) # don't call at nirvana since that doesn't correspond to a cell, now polys_clean_scc has the connected polygons, polys_clean is the rest
    polys_centers_clean_scc = polys_centers_clean[scc_inds[1:end-1]]
    reinjection_counts = reinjection_counts[scc_inds]

    size_change = initial_size - Psize

    if size_change == 0
        scc_info = "The graph is connected, no data was removed. The number of polygons is " * string(Psize - 1)
    else
        scc_info = "The graph is disconnected. Taking the strongest connected component killed " * string(size_change) * " components. The number of polygons is " * string(Psize - 1) * "."
    end

    ####################################################################################################################################
    """
    Collect polygons.
    """
    ##########################################################################################

    # Collect cells into n x 3 array of the form: cell # | x | y
    # Simultaneously, collect polygon centers into n x 2 array of the form: | x | y | 
    vcenters = zeros(length(polys_clean_scc), 2)
    vcells = zeros(sum(length(polys_clean_scc[i]) for i = 1:Psize-1), 3)
    v_cell_i = 1
    for n = 1:Psize-1
        for point in polys_clean_scc[n]
            vcells[v_cell_i, :] = [n, point[1], point[2]]
            v_cell_i = v_cell_i + 1
        end

        vcenters[n,:] = polys_centers_clean_scc[n]
    end

    # the set of disconnected polys
    if size_change == 0
        vcells_dis = []
    else
        vcells_dis = zeros(sum(length(polys_clean[i]) for i = 1:length(polys_clean)), 3)
        v_cell_i = 1
        for n = 1:length(polys_clean)
            for point in polys_clean[n]
                vcells_dis[v_cell_i, :] = [n, point[1], point[2]]
                v_cell_i = v_cell_i + 1
            end
        end
    end

    ####################################################################################################################################
    """
    Reinject trajectories which go from nirvana to the interior and normalize P.
    """
    ##########################################################################################

    if sto_type == "data"
        P_closed[Psize,:] = reinjection_counts
    elseif sto_type == "source"
        source_inds = getABinds(vcells, source_centers)
        P_closed[Psize, source_inds] .= 1
    end

    final_counts = [sum(P_closed[i,:]) for i = 1:Psize]

    # normalize P
    for i = 1:Psize
        if !(i == Psize && final_counts[i] == 0) # there are no transitions from nirvana if this is true. If it is, fine to leave zero row at end.
            P_closed[i,:] = P_closed[i,:]/final_counts[i]
        end
    end

    ####################################################################################################################################
    """
    Return results.
    """
    ##########################################################################################
    
    pi_closed = abs.(normalize(eigvecs(transpose(P_closed))[:,size(P_closed)[1]], 1))
    P_open = P_closed[1:end - 1, 1:end - 1]
    pi_open = pi_closed[1:end - 1]
    leaves = [(1.0 - P_open[i, i])*final_counts[i] for i = 1:Psize - 1]
        
    info = clean_info * scc_info

    returndict = begin Dict(
        "P_closed" => P_closed, 
        "P_open" => P_open,
        "pi_closed" => pi_closed,
        "pi_open" => pi_open,
        "counts" => final_counts,
        "leaves" => leaves,
        "polys" => vcells,
        "polys_centers" => vcenters,
        "polys_dis" => vcells_dis,
        "info" => info)
    end
    
    return returndict 
end

#refactor