var documenterSearchIndex = {"docs":
[{"location":"#UlamMethod.jl","page":"UlamMethod.jl","title":"UlamMethod.jl","text":"","category":"section"},{"location":"","page":"UlamMethod.jl","title":"UlamMethod.jl","text":"Documentation for UlamMethod.jl (test!)","category":"page"},{"location":"","page":"UlamMethod.jl","title":"UlamMethod.jl","text":"Modules = [UlamMethod]\nOrder = [:function, :type]","category":"page"},{"location":"#UlamMethod.binner_square-Tuple{UlamDomain}","page":"UlamMethod.jl","title":"UlamMethod.binner_square","text":"binner_square(domain)\n\nCover the computational domain in domain by a uniform grid of squares. Returns Vector{UlamPolygon{Float64}}.\n\n\n\n\n\n","category":"method"},{"location":"#UlamMethod.binner_voronoi-Tuple{UlamTrajectories, UlamDomain}","page":"UlamMethod.jl","title":"UlamMethod.binner_voronoi","text":"binner_voronoi(traj, domain; rseed)\n\nCover the domain by a Voronoi tesselation based on points generated by  kmeans clustering of the initial trajectory points. Returns Vector{UlamPolygon{Float64}}.\n\n\n\n\n\n","category":"method"},{"location":"#UlamMethod.inpoly-Tuple{Matrix{Float64}, PolyTable}","page":"UlamMethod.jl","title":"UlamMethod.inpoly","text":"inpoly(data::Matrix{Float64}, polys::PolyTable)\n\nDetermines which polygon of polys contains the data points in data. Returns an InpolyResult[@ref].\n\n\n\n\n\n","category":"method"},{"location":"#UlamMethod.inpoly-Tuple{UlamTrajectories, UlamDomain}","page":"UlamMethod.jl","title":"UlamMethod.inpoly","text":"inpoly(traj::UlamTrajectories, domain::UlamDomain)\n\nDetermines the indices of points that are inside the domain (i.e. not nirvana.)\n\n\n\n\n\n","category":"method"},{"location":"#UlamMethod.ulam_binner-Tuple{UlamTrajectories, UlamDomain}","page":"UlamMethod.jl","title":"UlamMethod.ulam_binner","text":"ulam_binner(traj, domain)\n\nSelects and executes the appropriate binning algorithm for the given trajectories and domain. Intersects the resulting polygons with the domain so that the returned polygons are clipped to it. Returns Vector{UlamPolygon{Float64}}.\n\n\n\n\n\n","category":"method"},{"location":"#UlamMethod.ulam_intersection-Tuple{UlamPolygon, UlamPolygon}","page":"UlamMethod.jl","title":"UlamMethod.ulam_intersection","text":"ulam_intersection(poly1, poly2)\n\nCompute the intersection of two UlamPolygons objects. Returns false if they do not intersect.\n\n\n\n\n\n","category":"method"},{"location":"#UlamMethod.ulam_intersects-Tuple{UlamPolygon, UlamPolygon}","page":"UlamMethod.jl","title":"UlamMethod.ulam_intersects","text":"ulam_intersects(poly1, poly2)\n\nReturn truw or false accoring to whether two UlamPolygons objects intersect. Can be faster than ulam_intersection if the shape of the intersection is not needed.\n\n\n\n\n\n","category":"method"},{"location":"#UlamMethod.ulam_method-Tuple{UlamTrajectories, UlamDomain}","page":"UlamMethod.jl","title":"UlamMethod.ulam_method","text":"ulam_method(traj, domain)\n\nThe high-level Ulam method; cover the domain by polygons and then construct the transition probability matrix.\n\n\n\n\n\n","category":"method"},{"location":"#UlamMethod.ulam_nirvana-Union{Tuple{T}, Tuple{UlamTrajectories, UlamDomain, Array{UlamPolygon{T}, 1}}} where T<:Real","page":"UlamMethod.jl","title":"UlamMethod.ulam_nirvana","text":"ulam_nirvana(traj, domain, polys)\n\nThis is the implementation of Ulam's method with a nirvana state. Returns an UlamResult object. Returns an UlamResult\n\n\n\n\n\n","category":"method"}]
}
