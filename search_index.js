var documenterSearchIndex = {"docs":
[{"location":"api/#API","page":"API","title":"API","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"These are the methods and types used by UlamMethod.jl.","category":"page"},{"location":"api/#Index","page":"API","title":"Index","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"Pages = [\"api.md\"]","category":"page"},{"location":"api/#Functions","page":"API","title":"Functions","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"Modules = [UlamMethod]\nOrder = [:function]","category":"page"},{"location":"api/#UlamMethod.binner_square-Tuple{UlamDomain}","page":"API","title":"UlamMethod.binner_square","text":"binner_square(domain)\n\nCover domain by a tight uniform grid of squares.\n\n\n\n\n\n","category":"method"},{"location":"api/#UlamMethod.binner_voronoi-Tuple{UlamTrajectories, UlamDomain}","page":"API","title":"UlamMethod.binner_voronoi","text":"binner_voronoi(traj, domain)\n\nCover the domain by a Voronoi tesselation based on points generated by  kmeans clustering of the initial trajectory points in traj.\n\n\n\n\n\n","category":"method"},{"location":"api/#UlamMethod.inpoly-Tuple{Matrix{Float64}, PolyTable}","page":"API","title":"UlamMethod.inpoly","text":"inpoly(data::Matrix{Float64}, polys::PolyTable)\n\nDetermine which polygon of polys contains the data points in data. Return an InpolyResult.\n\n\n\n\n\n","category":"method"},{"location":"api/#UlamMethod.inpoly-Tuple{UlamTrajectories, UlamDomain}","page":"API","title":"UlamMethod.inpoly","text":"inpoly(traj::UlamTrajectories, domain::UlamDomain)\n\nDetermine the indices of points that are inside domain.domain. Return an InpolyResult.\n\n\n\n\n\n","category":"method"},{"location":"api/#UlamMethod.ulam_binner-Tuple{UlamTrajectories, UlamDomain}","page":"API","title":"UlamMethod.ulam_binner","text":"ulam_binner(traj, domain)\n\nSelect and executes the appropriate binning algorithm for the given trajectories and domain.\n\nIntersect the resulting polygons with the domain so that the returned polygons are clipped to it.\n\n\n\n\n\n","category":"method"},{"location":"api/#UlamMethod.ulam_intersection-Tuple{UlamPolygon, UlamPolygon}","page":"API","title":"UlamMethod.ulam_intersection","text":"ulam_intersection(poly1, poly2)\n\nCompute the intersection of two UlamPolygon objects. Return false if they do not intersect.\n\n\n\n\n\n","category":"method"},{"location":"api/#UlamMethod.ulam_intersects-Tuple{UlamPolygon, UlamPolygon}","page":"API","title":"UlamMethod.ulam_intersects","text":"ulam_intersects(poly1, poly2)\n\nReturn true or false accoring to whether two UlamPolygon objects intersect. Can be faster than ulam_intersection if the shape of the intersection is not needed.\n\n\n\n\n\n","category":"method"},{"location":"api/#UlamMethod.ulam_method-Tuple{UlamTrajectories, UlamDomain}","page":"API","title":"UlamMethod.ulam_method","text":"ulam_method(traj, domain)\n\nRun the high-level Ulam method and return an UlamResult.\n\nArguments\n\ntraj: An UlamTrajectories; contains the trajectory data.\ndomain: An UlamDomain; contains the domain specification.\n\n\n\n\n\n","category":"method"},{"location":"api/#UlamMethod.ulam_nirvana-Union{Tuple{T}, Tuple{UlamTrajectories, UlamDomain, Array{UlamPolygon{T}, 1}}} where T<:Real","page":"API","title":"UlamMethod.ulam_nirvana","text":"ulam_nirvana(traj, domain, polys)\n\nExecute the core algorithm of Ulam's method with a nirvana state. Return an UlamResult.\n\n\n\n\n\n","category":"method"},{"location":"api/#UlamMethod.ulam_write-Tuple{String, UlamResult}","page":"API","title":"UlamMethod.ulam_write","text":"ulam_write(outfile, ulam_result; P_out)\n\nWrite ulam_result to the file outfile, which must be in the .h5 format. \n\nOptional Arguments\n\nP_out: If false, P_closed is not written to file since it can sometimes be very large. Default true.\n\n\n\n\n\n","category":"method"},{"location":"api/#Types","page":"API","title":"Types","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"Modules = [UlamMethod]\nOrder = [:type]","category":"page"},{"location":"api/#UlamMethod.InpolyResult","page":"API","title":"UlamMethod.InpolyResult","text":"InpolyResult{U}\n\nHold the result of inpoly.\n\n\n\n\n\n","category":"type"},{"location":"api/#UlamMethod.PolyTable","page":"API","title":"UlamMethod.PolyTable","text":"Polytable{T, U}\n\nTables, in matrix form of nodes and edges making up a set of polygons.\n\n\n\n\n\n","category":"type"},{"location":"api/#UlamMethod.PolyTable-Union{Tuple{Array{UlamPolygon{T}, 1}}, Tuple{T}} where T<:Real","page":"API","title":"UlamMethod.PolyTable","text":"Polytable(ulam_polys)\n\nConstruct a table of matrices and nodes of a vector of UlamPolygons.\n\n\n\n\n\n","category":"method"},{"location":"api/#UlamMethod.UlamDomain","page":"API","title":"UlamMethod.UlamDomain","text":"UlamDomain{S, T, U}\n\nThe domain defining the Ulam problem setup.\n\nFields\n\ncorners: The rectangular bounding box of the polygon in domain. Binning algorithms start by covering this rectangle.\ndomain: Data outside this are considered in nirvana.\npoly_type: The kind of polygons which cover domain.\npoly_number: The number of polygons which cover domain.\nstoc_type: The manner in which trajectories which leave domain are reinjected.\nstoc_source: The location in which to reinject data for the \"source\" reinjection algorithm.\n\n\n\n\n\n","category":"type"},{"location":"api/#UlamMethod.UlamDomain-Union{Tuple{UlamPolygon{T}}, Tuple{U}, Tuple{T}, Tuple{S}} where {S<:AbstractString, T<:Real, U<:Integer}","page":"API","title":"UlamMethod.UlamDomain","text":"UlamDomain(domain; [...])\n\nConstruct an UlamDomain defined by the UlamPolygon in domain.\n\nPoints outside domain will be considered in nirvana.\n\nOptional Arguments\n\npoly_type: One of \"sqr\", \"hex\", and \"vor\" for coverings by squares, hexagons or Voronoi tesselation. The default is squares.\npoly_number: The number of polygons requested. The default is 500 for squares/hexagons and 100 for Voronoi.\nstoc_type: Picks the stochasticization algorithm; one of \"data\" or \"source\". The default is data.\nstoc_source::UlamPolygon: Polygons in the covering that intersect stoc_source will have data re-injected uniformly through them in the source algorithm. \nstoc_source::Matrix: Polygons in the covering which contain the points in the stoc_source matrix will have data re-injected uniformly through them in the source algorithm. \nrseed: A seed for reproducing the random initialization of the kmeans algorithm in the Voronoi covering.\n\n\n\n\n\n","category":"method"},{"location":"api/#UlamMethod.UlamDomain-Union{Tuple{U}, Tuple{T}, Tuple{S}, NTuple{4, T}} where {S<:AbstractString, T<:Real, U<:Integer}","page":"API","title":"UlamMethod.UlamDomain","text":"UlamDomain(xmin, xmax, ymin, ymax; [...])\n\nConstruct an UlamDomain defined by the rectangle with bottom left corner (xmin, ymin) and top right corner (xmax, ymax).\n\nThis is equivalent to UlamDomain(domain; [...]) where the provided domain is a rectangular UlamPolygon.\n\n\n\n\n\n","category":"method"},{"location":"api/#UlamMethod.UlamInfo","page":"API","title":"UlamMethod.UlamInfo","text":"UlamInfo{S, U}\n\nA container for some of the results of ulam_method.\n\n\n\n\n\n","category":"type"},{"location":"api/#UlamMethod.UlamPolygon","page":"API","title":"UlamMethod.UlamPolygon","text":"UlamPolygon{T}\n\nA polygon defined by a series of nodes and carrying a center.\n\n\n\n\n\n","category":"type"},{"location":"api/#UlamMethod.UlamPolygon-Tuple{Matrix{<:Real}}","page":"API","title":"UlamMethod.UlamPolygon","text":"UlamPolygon(nodes; edges = nothing)\n\nConstruct an UlamPolygon based on a nodes matrix with n rows and 2 columns.\n\nThe edges of the polygon defined by nodes are connected in order.\n\nOptional Arguments\n\nedges: an n by 2 matrix which specifies the edge connections betwen nodes. Used if nodes are not already sorted.\n\n\n\n\n\n","category":"method"},{"location":"api/#UlamMethod.UlamResult","page":"API","title":"UlamMethod.UlamResult","text":"UlamResult{S, T, U}\n\nA container for the output of ulam_method.\n\n\n\n\n\n","category":"type"},{"location":"api/#UlamMethod.UlamResult-Union{Tuple{U}, Tuple{T}, Tuple{S}, Tuple{Matrix{T}, Array{UlamPolygon{T}, 1}, Array{UlamPolygon{T}, 1}, Vector{U}, UlamInfo{S, U}}} where {S<:AbstractString, T<:Real, U<:Integer}","page":"API","title":"UlamMethod.UlamResult","text":"UlamResult(P_closed, polys, polys_dis, counts, info)\n\nConstruct a container for the output of ulam_method.\n\nArguments\n\nP_closed: the transition probability matrix obtained by Ulam's method.\npolys: the vector of UlamPolygons constituting the covering of the domain.\npolys_dis: the vector of UlamPolygons which contained data but were in disconnected components.\ncounts: the number of observation data points in the ith polygon.\ninfo: additional info from UlamInfo.\n\n\n\n\n\n","category":"method"},{"location":"api/#UlamMethod.UlamTrajectories","page":"API","title":"UlamMethod.UlamTrajectories","text":"UlamTrajectories{T}\n\nA container for trajectory data.\n\n\n\n\n\n","category":"type"},{"location":"api/#UlamMethod.UlamTrajectories-Tuple{String}","page":"API","title":"UlamMethod.UlamTrajectories","text":"UlamTrajectories(infile, [...])\n\nConstruct a container for trajectory data, loading it from infile, which should be a .h5 or .mat file.\n\nOptional Arguments\n\nThese are used if trajectory data in infile are named something other than \"x0\", \"y0\", \"xT\", \"yT\".\n\nx0_alias\ny0_alias\nxT_alias\nyT_alias\n\n\n\n\n\n","category":"method"},{"location":"api/#UlamMethod.UlamTrajectories-Tuple{}","page":"API","title":"UlamMethod.UlamTrajectories","text":"UlamTrajectories(;x0, y0, xT, yT)\n\nConstruct a container for trajectory data.\n\n\n\n\n\n","category":"method"},{"location":"api/#Constants","page":"API","title":"Constants","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"Modules = [UlamMethod]\nOrder = [:constant]","category":"page"},{"location":"api/#UlamMethod.global_poly_number_default","page":"API","title":"UlamMethod.global_poly_number_default","text":"global_poly_number_default\n\nThe default number of polygons for each of the types of coverings.\n\n\"sqr\": 500\n\"hex\": 500\n\"vor\": 100\n\n\n\n\n\n","category":"constant"},{"location":"api/#UlamMethod.global_poly_types","page":"API","title":"UlamMethod.global_poly_types","text":"global_poly_types\n\nOne of \"sqr\", \"hex\" and \"vor\" according to the coverings in UlamDomain.\n\n\n\n\n\n","category":"constant"},{"location":"api/#UlamMethod.global_rseed_default","page":"API","title":"UlamMethod.global_rseed_default","text":"global_rseed_default\n\nThe default seed for RNG reproducibility. Default 123.\n\n\n\n\n\n","category":"constant"},{"location":"api/#UlamMethod.global_stoc_types","page":"API","title":"UlamMethod.global_stoc_types","text":"global_stoc_types\n\nOne of \"data\" or \"source\" according to the stochasticization options in UlamDomain.\n\n\n\n\n\n","category":"constant"},{"location":"api/#UlamMethod.global_traj_file_types","page":"API","title":"UlamMethod.global_traj_file_types","text":"global_traj_file_types\n\nOne of \"mat\", \"h5\" according to the trajectory input file in UlamTrajectories\n\n\n\n\n\n","category":"constant"},{"location":"api/","page":"API","title":"API","text":"","category":"page"},{"location":"theory/#Theory-and-Implementation","page":"Theory and Implementation","title":"Theory and Implementation","text":"","category":"section"},{"location":"theory/","page":"Theory and Implementation","title":"Theory and Implementation","text":"UNDER CONSTRUCTION.","category":"page"},{"location":"theory/#Calculation-Hierarchy","page":"Theory and Implementation","title":"Calculation Hierarchy","text":"","category":"section"},{"location":"theory/","page":"Theory and Implementation","title":"Theory and Implementation","text":"The user provides UlamTrajectories via file or manual input. This type contains x0, y0, xT, yT.\nThe user provides an UlamDomain which contains the domain, an UlamPolygon defining the boundary of nirvana. Also provided is the type and number of polygons requested, and the type and location of the reinjection algorithm.\nThe highest level function is ulam_method which calls the various steps in the calculation.\nFirst, ulam_binner is used. This selects and runs the appropriate binning algorithm based on the user's choice. Every binning algorithm acts on the bounding box of domain. The output of every binning algorithm is vector of UlamPolygons which represent a covering of the bounding box.\nAfterwards, the domain is intersected with the binning result, again giving a vector of UlamPolygon but such that every data point is now inside domain. The intersection is such that if the result of an intersection is a multipolygon, the polygon with the largest area is taken.\nThe result of ulam_binner \"actual\" Ulam's method calculation, which is contained in ulam_nirvana. This function takes the trajectories and polygons and constructs the transition matrix defining the Markov chain on states representing the polygons. Polygons which contain no (x0, y0) points are discarded. The largest strongly connected component of the transition matrix is extracted. Both of these \"cleaning\" operations have the effect of deleting trajectories which start or end in any of the removed polygons.","category":"page"},{"location":"advanced/#Advanced-Usage","page":"Advanced Usage","title":"Advanced Usage","text":"","category":"section"},{"location":"advanced/#Loading-trajectories","page":"Advanced Usage","title":"Loading trajectories","text":"","category":"section"},{"location":"advanced/","page":"Advanced Usage","title":"Advanced Usage","text":"For larger systems, trajectory data may be loaded from a file. The core functionality is provided by ","category":"page"},{"location":"advanced/","page":"Advanced Usage","title":"Advanced Usage","text":"UlamTrajectories(infile; x0_alias, y0_alias, xT_alias, yT_alias)","category":"page"},{"location":"advanced/","page":"Advanced Usage","title":"Advanced Usage","text":"This file should be in the .mat or .h5 format with the keys \"x0\", \"y0\", \"xT\" and \"xT\" in the root of the file.  Aliases for the trajectory data may be provided. For example, if we had a file in which the \"x0\" data were instead called \"xtraj0\", we could use traj = UlamTrajectories(infile, x0_alias = \"xtraj0\") and so on for other keys.","category":"page"},{"location":"advanced/#Defining-domains","page":"Advanced Usage","title":"Defining domains","text":"","category":"section"},{"location":"advanced/","page":"Advanced Usage","title":"Advanced Usage","text":"The core functionality is provided by ","category":"page"},{"location":"advanced/","page":"Advanced Usage","title":"Advanced Usage","text":"UlamDomain(domain; poly_type, poly_number, stoc_type, stoc_polygon, rseed)","category":"page"},{"location":"advanced/","page":"Advanced Usage","title":"Advanced Usage","text":"The field domain should be an UlamPolygon. All data outside domain will be considered to be in nirvana.[5] [6] A convenience method is provided for square domains:","category":"page"},{"location":"advanced/","page":"Advanced Usage","title":"Advanced Usage","text":"UlamDomain(xmin, xmax, ymion, ymax; poly_type, poly_number, stoc_type, stoc_polygon, rseed)","category":"page"},{"location":"advanced/","page":"Advanced Usage","title":"Advanced Usage","text":"This automatically constructs a rectangular domain with bottom left corner (xmin, ymin) and upper right corner (xmax, ymax). In either case, the rectangular bounding box of the domain is automatically computed. Refer to src/earth-polygons.jl for some sample domains.","category":"page"},{"location":"advanced/#binning","page":"Advanced Usage","title":"Binning Algorithms","text":"","category":"section"},{"location":"advanced/","page":"Advanced Usage","title":"Advanced Usage","text":"These algorithms control how the computational rectangle (that is, the bounding box of the domain) is covered in polygons. The named argument poly_type to UlamDomain can take one of three values:","category":"page"},{"location":"advanced/","page":"Advanced Usage","title":"Advanced Usage","text":"\"sqr\": A covering of the domain by a regular grid of squares (default.) This algorithm is suitable for trajectory data with few \"holes\" in the observations.\n\"hex\": A covering of the domain by a regular grid of hexagons. Hexagons may have some advantages over squares [1].\n\"vor\": A covering of the domain by a Voronoi tesselation generated by cluster centers obtained from k-means [2]. This algorithm is more appropriate for sparse data. Pass rseed for reproducible results; the kmeans initialization contains some randomness.","category":"page"},{"location":"advanced/","page":"Advanced Usage","title":"Advanced Usage","text":"The named argument poly_number selects the number of polygons in each case.","category":"page"},{"location":"advanced/#stoc","page":"Advanced Usage","title":"Stochasticization Algorithms","text":"","category":"section"},{"location":"advanced/","page":"Advanced Usage","title":"Advanced Usage","text":"These algorithms control how reinjection counts (trajectories pointing from nirvana to the interior) are distributed. The named argumet stoc_type to UlamDomain can take one of two values:","category":"page"},{"location":"advanced/","page":"Advanced Usage","title":"Advanced Usage","text":"\"data\": Reinjection occurs according to which boxes trajectories actually enter (this is the default algorithm.)\n\"source\": Reinjection occurs uniformly in polygons specified by stoc_source. If stoc_source is provided, \"source\" is automatically selected. If \"source\" is selected but no stoc_source is provided, then stoc_source is set equal to the domain. This is equivalent to reinjecting data uniformly across all boxes.\nIf stoc_source is entered as an N times 2 matrix of numbers, then reinjection occurs at the set of polygons which contain at least one point from stoc_source.\nIf stoc_source is entered as an UlamPolygon, then reinjection occurs at the set of polygons which intersect stoc_source.","category":"page"},{"location":"advanced/#Using-the-results","page":"Advanced Usage","title":"Using the results","text":"","category":"section"},{"location":"advanced/","page":"Advanced Usage","title":"Advanced Usage","text":"The core functionality is provided by ","category":"page"},{"location":"advanced/","page":"Advanced Usage","title":"Advanced Usage","text":"ulam_method(traj, domain)","category":"page"},{"location":"advanced/","page":"Advanced Usage","title":"Advanced Usage","text":"The output of ulam_method is an UlamResult. This contains the main objects calculated by Ulam's Method.","category":"page"},{"location":"advanced/","page":"Advanced Usage","title":"Advanced Usage","text":"polys: A vector of UlamPolygons which define the covering. Use PolyTable to access the nodes in a readable format.\npolys_dis: A vector of UlamPolygons which contained data but were disconnected when the strongest connected component was calculated.\nP_closed: The full transition matrix. The last row and column correspond to nirvana.\npi_closed: The largest left eigenvector of P_closed.","category":"page"},{"location":"advanced/","page":"Advanced Usage","title":"Advanced Usage","text":"For convenience P_open and pi_open are also provided, which are identical to P_closed and pi_closed with the nirvana entries removed.","category":"page"},{"location":"advanced/#Writing-the-results","page":"Advanced Usage","title":"Writing the results","text":"","category":"section"},{"location":"advanced/","page":"Advanced Usage","title":"Advanced Usage","text":"A custom write method is provided","category":"page"},{"location":"advanced/","page":"Advanced Usage","title":"Advanced Usage","text":"ulam_write(outfile, ulam_result; P_out)","category":"page"},{"location":"advanced/","page":"Advanced Usage","title":"Advanced Usage","text":"This will write an UlamResult to the file specified by outfile. Note that outfile must be of the form \"fname.h5\". Optionally pass P_out = false to avoid writing the P_closed matrix since it can be very large.","category":"page"},{"location":"advanced/","page":"Advanced Usage","title":"Advanced Usage","text":"The polygons will be output in an N times 3 matrix such that the first two columns are the (x y) coordinates of a polygon vertex and the third column is the index of the polygon that vertex belongs to. The vertices are sorted. ","category":"page"},{"location":"advanced/#Full-workflow-example","page":"Advanced Usage","title":"Full workflow example","text":"","category":"section"},{"location":"advanced/","page":"Advanced Usage","title":"Advanced Usage","text":"For this example, the file test/x0x5-NA-undrogued.h5 contains trajectory data from undrogued drifters in the North Atlantic obtained from the NOAA GDP [3] [4].","category":"page"},{"location":"advanced/","page":"Advanced Usage","title":"Advanced Usage","text":"infile = \"x0x5-NA-undrogued.h5\"     # place this file in your working directory, or define a path to it\ntraj = UlamTrajectories(infile)","category":"page"},{"location":"advanced/","page":"Advanced Usage","title":"Advanced Usage","text":"Next we define our domain. We'll use North_Atlantic_clipped_verts here. For the binning, we'll use the default square covering with 760 polygons. We'll also use the default \"data\" stochasticization algorithm.","category":"page"},{"location":"advanced/","page":"Advanced Usage","title":"Advanced Usage","text":"NA = UlamPolygon(North_Atlantic_clipped_verts)\npoly_type = \"sqr\"\npoly_number = 760\n\ndomain = UlamDomain(NA, poly_type = poly_type, poly_number = poly_number)","category":"page"},{"location":"advanced/","page":"Advanced Usage","title":"Advanced Usage","text":"The final step is to apply Ulam's method.","category":"page"},{"location":"advanced/","page":"Advanced Usage","title":"Advanced Usage","text":"ulam = ulam_method(traj, domain)","category":"page"},{"location":"advanced/","page":"Advanced Usage","title":"Advanced Usage","text":"From here, a .h5 file can be created with ulam_write(\"my_ulam_results.h5\", ulam) or ulam can be used elsewhere.","category":"page"},{"location":"advanced/#References","page":"Advanced Usage","title":"References","text":"","category":"section"},{"location":"advanced/","page":"Advanced Usage","title":"Advanced Usage","text":"[1]: https://www.uber.com/blog/h3/","category":"page"},{"location":"advanced/","page":"Advanced Usage","title":"Advanced Usage","text":"[2]: https://juliastats.org/Clustering.jl/stable/kmeans.html#K-means","category":"page"},{"location":"advanced/","page":"Advanced Usage","title":"Advanced Usage","text":"[3]: https://www.aoml.noaa.gov/phod/gdp/data.php","category":"page"},{"location":"advanced/","page":"Advanced Usage","title":"Advanced Usage","text":"[4]: Lumpkin, Rick, and Mayra Pazos. \"Measuring surface currents with Surface Velocity Program drifters: the instrument, its data, and some recent results.\" Lagrangian analysis and prediction of coastal and ocean dynamics 39 (2007): 67.","category":"page"},{"location":"advanced/","page":"Advanced Usage","title":"Advanced Usage","text":"[5]: Miron, Philippe, et al. \"Transition paths of marine debris and the stability of the garbage patches.\" Chaos: An Interdisciplinary Journal of Nonlinear Science 31.3 (2021): 033101.","category":"page"},{"location":"advanced/","page":"Advanced Usage","title":"Advanced Usage","text":"[6]: In brief, nirvana is an extra state appended to an open system to close it; trajectories which point from inside the domain to the outisde of the domain transition to this nirvana state. Trajectories which point from outside the domain to the inside are transitions \"from\" nirvana - how exactly these data are reinjected is controlled by the stochasticization algorithms.","category":"page"},{"location":"#UlamMethod.jl","page":"Home","title":"UlamMethod.jl","text":"","category":"section"},{"location":"#Introduction","page":"Home","title":"Introduction","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"This package is an implementation of Ulam's method [1] [2] (see also Galerkin projection [3]) for the discretization of a stochastic operator using pure Julia. Given a set of two-dimensional, one-step trajectories ","category":"page"},{"location":"","page":"Home","title":"Home","text":"(x_0 1 y_0 1) to  (x_T 1 y_T 1) (x_0 2 y_0 2) to  (x_T 2 y_T 2) dots","category":"page"},{"location":"","page":"Home","title":"Home","text":"defined in a domain, the essential goal of Ulam's method is to partition the domain into a series of non-intersecting boxes and construct a transition probability matrix P on these boxes. In UlamMethod.jl, this is accomplished in two main steps","category":"page"},{"location":"","page":"Home","title":"Home","text":"The user provides a polygon containing the data, and covering of the the domain is generated by polygons according to one of several different binning algorithms.\nThe number of trajectories beginning in polygon i and ending in polygon j is used to create the entry P_i j of P such that the edges of the domain are handled by a stochasticization algorithm.","category":"page"},{"location":"","page":"Home","title":"Home","text":"The polygons which form the covering and the transition probability matrix are the main outputs.","category":"page"},{"location":"#Installation","page":"Home","title":"Installation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"In the Julia REPL, run the following code and follow the prompts:","category":"page"},{"location":"","page":"Home","title":"Home","text":"using Pkg\nPkg.add(\"UlamMethod\")","category":"page"},{"location":"","page":"Home","title":"Home","text":"Make the functions in this package available to use in your code by including the following line:","category":"page"},{"location":"","page":"Home","title":"Home","text":"using UlamMethod","category":"page"},{"location":"#Quickstart","page":"Home","title":"Quickstart","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The core functionality is provided by ","category":"page"},{"location":"","page":"Home","title":"Home","text":"ulam_method(traj::UlamTrajectories, domain::UlamDomain)","category":"page"},{"location":"","page":"Home","title":"Home","text":"where traj contains information about trajectories and domain contains information about the domain and covering. Here are 10000 random trajectories in the domain 0 10^2","category":"page"},{"location":"","page":"Home","title":"Home","text":"n_data = 10000\nx0, y0, xT, yT = 10*rand(n_data), 10*rand(n_data), 10*rand(n_data), 10*rand(n_data)\ntraj = UlamTrajectories(x0 = x0, y0 = y0, xT = xT, yT = yT)","category":"page"},{"location":"","page":"Home","title":"Home","text":"We will take our domain to be the rectangular subset 3 5 times 4 8 and generate a covering with 40 squares.","category":"page"},{"location":"","page":"Home","title":"Home","text":"xmin, xmax, ymin, ymax = 3, 5, 4, 8\ndomain = UlamDomain(xmin, xmax, ymin, ymax, poly_type = \"sqr\", poly_number = 40)","category":"page"},{"location":"","page":"Home","title":"Home","text":"Run Ulam's method.","category":"page"},{"location":"","page":"Home","title":"Home","text":"ulam = ulam_method(traj, domain)    # the main calculation\n\nulam.P_closed                       # the transition matrix\npt = PolyTable(ulam.polys)          # PolyTable makes a simple list of nodes and edges\npt.nodes                            # |x|y| table of polygon nodes\npt.edges[:,3]                       # the index of the polygon that the i'th node belongs to","category":"page"},{"location":"#References","page":"Home","title":"References","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"[1]: Ulam, Stanislaw M. A collection of mathematical problems. No. 8. Interscience Publishers, 1960.","category":"page"},{"location":"","page":"Home","title":"Home","text":"[2]: Li, Tien-Yien. \"Finite approximation for the Frobenius-Perron operator. A solution to Ulam's conjecture.\" Journal of Approximation theory 17.2 (1976): 177-186.","category":"page"},{"location":"","page":"Home","title":"Home","text":"[3]: Reddy, Junuthula Narasimha. Introduction to the finite element method. McGraw-Hill Education, 2019.","category":"page"}]
}
