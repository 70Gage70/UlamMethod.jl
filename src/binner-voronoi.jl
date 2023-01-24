"""
Bin the data based on a kmeans + Voronoi tesselation.

The main steps are:

1. Cluster the observations using k-means, handled by the functions kmeans in the Julia package Clustering
2. Construct the polygons defining the boxes using a Voronoi tesselation, handled by the function voronoicells() in the Julia package VaronoiCells
3. Truncate edge polygons to ther intersection with the concave hull of the data (obs + traj).
4. Assign trajectory data xT to tesselation with inpoly2 from PolygonInbounds

Written by Gage Bonner November 2022
"""
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

using Random
using Clustering: kmeans
using VoronoiCells: voronoicells, Rectangle
using GeometryBasics: Point2
using LibGEOS
using GeoInterface
using PolygonInbounds
using LazySets: convex_hull

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

"""
Find indices of obs/traj data according to which are inside vs. outside nirvana. Return dictionary.
"""

function nirv_inds(x0, y0, xT, yT, corners)
    xmin, xmax, ymin, ymax = corners
    in0 = Int64[]
    inT = Int64[]
    nirv0 = Int64[]
    nirvT = Int64[]

    for i = 1:lastindex(x0)
        if (xmin <= x0[i] <= xmax) && (ymin <= y0[i] <= ymax)
            push!(in0, i)
        else
            push!(nirv0, i)
        end

        if (xmin <= xT[i] <= xmax) && (ymin <= yT[i] <= ymax)
            push!(inT, i)
        else
            push!(nirvT, i)
        end
    end
    
    return Dict("in0" => in0, "inT" => inT, "nirv0" => nirv0, "nirvT" => nirvT)
end

"""
Generates tesselation in bounding box defined by the variable corners with points "centers" (assumes centers is directly from kmeans).
The voronoicells function needs a rectangle defined by the function Rectangle and points constructed with Point2. These
are inside VoronoiCells.
"""

function voronoi(centers, corners)
    voronoi_rect = Rectangle(Point2(corners[1], corners[3]), Point2(corners[2], corners[4])) # xmin, ymin, xmax, ymax. "Rectangle" is from VoronoiCells
    voronoi_points = [Point2(centers[:, i]) for i = 1:length(centers[1, :])] # Point2 is from GeometryBasics
    voronoi_tess = voronoicells(voronoi_points, voronoi_rect).Cells
    voronoi_vec = [] # want to convert from Point2 etc. to just normal Floats 
    for i = 1:length(voronoi_tess)
        vc = [[voronoi_tess[i][j].data[1], voronoi_tess[i][j].data[2]] for j = 1:length(voronoi_tess[i])]
        push!(vc, [voronoi_tess[i][1].data[1], voronoi_tess[i][1].data[2]]) # close polygon by putting first point back
        push!(voronoi_vec, vc)
    end

    return voronoi_vec
end

"""
Converts from a polygon defined by a list of coordinates to Well-Known-Text (WKT) which is suitable for LibGEOS
"""

function coords_to_LibGEOS_polygon(coords)
    wkt = "POLYGON(("

    for i = 1:length(coords) - 1
        this_point = string(coords[i][1]) * " " * string(coords[i][2]) * ","
        wkt = wkt * this_point
    end

    wkt = wkt * string(coords[end][1]) * " " * string(coords[end][2]) * "))"

    return LibGEOS.readgeom(wkt)

end

"""
Computes the intersection of the Voronoi tesselation above with the convex hull of the data (requires the hull as a set of coordinates).
The output of voronoi_hull is a vector of length n_clusters, where the i'th entry of the vector is a vector of [x, y] vertices defining the i'th cell
"""

function voronoi_hull(centers, corners, hull)
    vt = voronoi(centers, corners)
    hull_poly = coords_to_LibGEOS_polygon(hull)

    new_vt = []

    for i = 1:length(vt)
        v_poly = coords_to_LibGEOS_polygon(vt[i])
        new_v = GeoInterface.coordinates(LibGEOS.intersection(hull_poly, v_poly))[1] # the actual intersection; LibGEOS "does" it and GeoInterface gets the coords
        push!(new_vt, new_v)
    end

    return new_vt
end

function voronoi_binner(x0, y0, xT, yT, n_clusters, corners; rseed = 123)
    nirv_dict = nirv_inds(x0, y0, xT, yT, corners)
    ####################################################################################################################################
    """
    Cluster observations (that are inside the computational domain) with kmeans.
    The kmeans algorithm takes a 2 x n matrix of data (so, want x0 to be the first row and y0 to be the second row)
    if res = kmeans(data, n_clusters), then res.clusters gives a 2 x n_clusters matrix of cluster centers and 
    res.assignments gives a column vector of assignments of data points to those clusters
    """
    println("Running kmeans.")
    data0 = transpose([x0[nirv_dict["in0"]] ;; y0[nirv_dict["in0"]]])
    Random.seed!(rseed) # random seed for kmeans algo; makes results reproducible
    res_kmeans = @time kmeans(data0, n_clusters)
    cluster_centers = res_kmeans.centers
    println("Done kmeans.")
    ####################################################################################################################################
    """
    Tile plane based on above clusters with Voronoi cells
    The Voronoi cells algorithm takes in an n-vector of Point2(x, y) as well as a rectangle that contains the points specified by Rectangle
    In our case, the rectangle should be defined by the variable corners since the kmeans were only clustered on data that's inside the box defined by corners
    """
    # println("Running Voronoi tesselation.")

    # First compute the convex hull of the data (observations + trajectories)
    convex_data = [[x0[i], y0[i]] for i in nirv_dict["in0"]] 
    convex_data = [convex_data ; [[xT[i], yT[i]] for i in nirv_dict["inT"]]]
    hull  = convex_hull(convex_data)
    push!(hull, hull[1]) # close polygon

    # println("Voronoi complete.")

    polys = voronoi_hull(cluster_centers, corners, hull)
    polys_centers = []

    for i = 1:n_clusters
        push!(polys_centers, cluster_centers[:,i])
    end

    return polys, polys_centers
end
