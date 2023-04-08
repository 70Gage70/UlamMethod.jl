using Random
using Clustering: kmeans
using VoronoiCells: voronoicells, Rectangle
using GeometryBasics: Point2
using Dates

"""
    binner_voronoi(traj, domain; rseed)

Cover the domain by a Voronoi tesselation based on points generated by 
kmeans clustering of the initial trajectory points.
Returns Vector{UlamPolygon{Float64}}.
"""
function binner_voronoi(traj::UlamTrajectories, domain::UlamDomain)
    # We only want to bin data inside the domain
    inds_in = inpoly(traj, domain).inds
    inds_in = [i == 1 ? true : false for i in inds_in]
    data0 = transpose([traj.x0[inds_in] ;; traj.y0[inds_in]])  # kmeans takes a 2 x n matrix; have to transpose

    # cluster the data
    tstart = Dates.format(now(), "HH:MM")
    @info "Starting kmeans with $(domain.poly_number) clusters at $(tstart)."
    Random.seed!(domain.rseed) # random seed for kmeans algo; makes results reproducible
    centers = @time kmeans(data0, domain.poly_number).centers
    @info "Done kmeans."

    # create tiling
    # The Voronoi cells algorithm takes in an n-vector of Point2(x, y) 
    # as well as a rectangle that contains the points specified by Rectangle
    xmin, xmax, ymin, ymax = domain.corners
    voronoi_rect = Rectangle(Point2(xmin, ymin), Point2(xmax, ymax)) # "Rectangle" is from VoronoiCells, Point2 is from GeometryBasics
    voronoi_points = [Point2(centers[:, i]) for i = 1:length(centers[1, :])]
    voronoi_tess = voronoicells(voronoi_points, voronoi_rect).Cells # vector of Vector{Vector}
    polys = [UlamPolygon(Matrix(reduce(hcat, vt)')) for vt in voronoi_tess]

    return polys
end
