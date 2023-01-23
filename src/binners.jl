"""
Bins 2D data according to binner algorithms.

Written by Gage Bonner December 2022.
"""


########################################################################################

include("binner-voronoi.jl")
include("binner-hexbin.jl")
include("binner-square.jl")

########################################################################################

using PolygonInbounds:inpoly2

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

function binner(x0, y0, type, n_polys; rseed = 123)
    if type == "vor"
        polys, polys_centers = voronoi_binner(x0, y0, n_polys, corners, rseed = rseed) 
    elseif type == "hex"
        polys, polys_centers = hexbin_binner(n_polys, corners) 
    elseif type == "reg"        
        polys, polys_centers = square_binner(n_polys, corners) 
    end

    n_polys_actual = length(polys)
    data0 = [x0 ;; y0]
    res = inpoly_preprocess(polys)
    nodes_inpoly = res["nodes"]
    edges_inpoly = res["edges"]

    inds0 = zeros(Int64, length(x0))
    contains_data = [] # incides of polygons that actually contain data

    res_poly0 = inpoly2(data0, nodes_inpoly, edges_inpoly)

    new_poly_count = 1 # keeping track of how many boxes actually contain data

    for i = 1:n_polys_actual
        res0i = res_poly0[:,1,i]
        if 1 in res0i # if there's at least one data point in this box
            push!(contains_data, i)
            these_inds0 = range(1, length(x0))[res0i]
            inds0[these_inds0] .= new_poly_count
            new_poly_count = new_poly_count + 1
        end
    end

    returnDict = Dict(
        "polys" => polys,
        "polys_centers" => polys_centers,
        "inds" => inds0,
        "contains_data" => contains_data
    )

    return returnDict
end