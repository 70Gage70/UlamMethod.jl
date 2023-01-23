"""
A generic Ulam method.
Needs: x0, y0, xT, yT, polys
Needs: polys: a vector whose i'th entry is a vector of [x, y] coordinates of that polygon WITH CLOSING NODE

Written by Gage Bonner November 2022
"""
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

using PolygonInbounds:inpoly2
using Graphs
using LinearAlgebra

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

"""
Generates a list of nodes and edges suitable for use in inpoly2()
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
The Ulam method. 

polys should be a vector of vectors; that is, poly[i] is a vector whose j'th entry is the jth vertex of that polygon. Must be closed!
polys_centers is also a vector of vectors; polys_centers[i] is equal to [x_center, y_center] 
"""

function ulam_nirvana(x0, y0, xT, yT, polys, polys_centers)
    ##########################################################################################
    """
    Assign indices to obs/traj data based on given polys.
    Simultaneously purge polygons with no observational data in them.
    """
    ##########################################################################################
    # println("Assigning indices.")

    npolys = length(polys)
    data0 = [x0 ;; y0]
    dataT = [xT ;; yT]
    res = inpoly_preprocess(polys)
    nodes_inpoly = res["nodes"]
    edges_inpoly = res["edges"]

    inds0 = zeros(Int64, length(x0))
    indsT = zeros(Int64, length(x0))
    contains_data = [] # incides of polygons that actually contain data

    res_poly0 = inpoly2(data0, nodes_inpoly, edges_inpoly)
    res_polyT = inpoly2(dataT, nodes_inpoly, edges_inpoly)

    new_poly_count = 1 # keeping track of how many boxes actually contain data

    for i = 1:npolys
        these_indsT = range(1, length(x0))[res_polyT[:,1,i]]
        indsT[these_indsT] .= new_poly_count

        res0i = res_poly0[:,1,i]
        if 1 in res0i # if there's at least one data point in this box
            push!(contains_data, i)
            these_inds0 = range(1, length(x0))[res0i]
            inds0[these_inds0] .= new_poly_count
            new_poly_count = new_poly_count + 1
        end
    end

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

    for i = 1:length(inds0)
        ind0 = inds0[i]
        indT = indsT[i]

        if ind0 != 0 && indT != 0
            P_closed[ind0, indT] = P_closed[ind0, indT] + 1
        elseif ind0 == 0 && indT != 0
            P_closed[Psize, indT] = P_closed[Psize, indT] + 1
        elseif ind0 != 0 && indT == 0
            P_closed[ind0, Psize] = P_closed[ind0, Psize] + 1
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
    polys_clean_scc = splice!(polys_clean, scc_inds[1:end-1]) # don't call at nirvana since that doesn't correspond to a cell, now vt_scc has the connected polygons, vt is the rest
    polys_centers_clean_scc = polys_centers_clean[scc_inds[1:end-1]]

    size_change = initial_size - Psize

    if size_change == 0
        scc_info = "The graph is connected, no data was removed. The number of polygons is " * string(Psize - 1)
    else
        scc_info = "The graph is disconnected. Taking the strongest connected component killed " * string(size_change) * " components. The number of polygons is " * string(Psize - 1) * "."
    end

    ####################################################################################################################################
    """
    Normalize P and return results.
    """
    ##########################################################################################
    # println("Cleaning and returning.")

    final_counts = [sum(P_closed[i,:]) for i = 1:Psize]

    # normalize P
    for i = 1:Psize
        if !(i == Psize && final_counts[i] == 0) # there are no transitions to nirvana if this is true. If it is, fine to leave zero row at end.
            P_closed[i,:] = P_closed[i,:]/final_counts[i]
        end
    end
    
    pi_closed = abs.(normalize(eigvecs(transpose(P_closed))[:,size(P_closed)[1]], 1))
    P_open = P_closed[1:end - 1, 1:end - 1]
    pi_open = pi_closed[1:end - 1]
    leaves = [(1.0 - P_open[i, i])*final_counts[i] for i = 1:Psize - 1]

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
        
    info = clean_info * scc_info

    returndict = begin Dict(
        "P_closed" => P_closed, 
        "P_open" => P_open,
        "pi_closed" => pi_closed,
        "pi_open" => pi_open,
        "counts" => final_counts,
        "leaves" => leaves,
        "polys" => vcells,
        "polys_raw" => polys_clean_scc,
        "polys_centers" => vcenters,
        "polys_dis" => vcells_dis,
        "info" => info)
    end
    
    return returndict 
end


function getABinds(ulam_polys_raw, A_centers, B_centers)
    res = inpoly_preprocess(ulam_polys_raw)
    nodes_inpoly = res["nodes"]
    edges_inpoly = res["edges"]

    indsA = zeros(Int64, length(ulam_polys_raw))
    indsB = zeros(Int64, length(ulam_polys_raw))

    res_polyA = inpoly2(A_centers, nodes_inpoly, edges_inpoly)
    res_polyB = inpoly2(B_centers, nodes_inpoly, edges_inpoly)

    for i = 1:length(ulam_polys_raw)
        if (1 in res_polyA[:,1,i]) # the i'th cell contains at least one point from A, therefore this cell belongs to A
            indsA[i] = 1
        end

        if (1 in res_polyB[:,1,i])
            indsB[i] = 1
        end

    end

    indsA = findall(!iszero, indsA)
    indsB = findall(!iszero, indsB)
    
    return Dict("indsA" => indsA, "indsB" => indsB)
end