"""
Bins data for the purposes of seeing its density. Uses squares.

Written by Gage Bonner December 2022
"""

using PolygonInbounds:inpoly2

include("binner-square.jl")
include("helpers.jl")

########################################################################################

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

function data_density(fname, nboxes)
    data_T = MAT.matopen("$fname.mat")
    x0, xT, y0, yT = read(data_T, "x0", "xT", "y0", "yT")
    x0 = vec(x0)
    y0 = vec(y0)
    xT = vec(xT)
    yT = vec(yT)
    close(data_T)

    corners = [-100, 15, -9, 39]

    polys, polys_centers = square_binner(nboxes, corners) 


    npolys = length(polys)
    data0 = [x0 ;; y0]
    res = inpoly_preprocess(polys)
    nodes_inpoly = res["nodes"]
    edges_inpoly = res["edges"]

    inds0 = zeros(Int64, length(x0))
    contains_data = [] # incides of polygons that actually contain data

    res_poly0 = inpoly2(data0, nodes_inpoly, edges_inpoly)

    new_poly_count = 1 # keeping track of how many boxes actually contain data

    for i = 1:npolys
        res0i = res_poly0[:,1,i]
        if 1 in res0i # if there's at least one data point in this box
            push!(contains_data, i)
            these_inds0 = range(1, length(x0))[res0i]
            inds0[these_inds0] .= new_poly_count
            new_poly_count = new_poly_count + 1
        end
    end

    polys_clean = polys[contains_data] # only keep the polygons that actually have data in them
    npolys = length(contains_data)

    counts = zeros(Int64, npolys)

    for ind in inds0
        if ind !=0 
            counts[ind] += 1
        end
    end

    # Collect cells into n x 3 array of the form: cell # | x | y
    # Simultaneously, collect polygon centers into n x 2 array of the form: | x | y | 
    vcells = zeros(sum(length(polys_clean[i]) for i = 1:npolys), 3)
    v_cell_i = 1
    for n = 1:npolys
        for point in polys_clean[n]
            vcells[v_cell_i, :] = [n, point[1], point[2]]
            v_cell_i = v_cell_i + 1
        end
    end

    csv("data_polys_" * string(nboxes), vcells);
    csv("data_counts_" * string(nboxes), counts);

    return "Done."
end