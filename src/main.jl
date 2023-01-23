using MAT
using HDF5

include("binner-voronoi.jl")
include("binner-hexbin.jl")
include("binner-square.jl")

include("ulam-nirvana.jl")
include("tpt-infinite.jl")

include("helpers.jl")
include("earth-polygons.jl")


########################################################################################
# SINGLE BOX A AND B
A_centers = [
    -18.0 17.0;
    ]
# B_centers = [
#     -95.0 27.0;
#     ]

# MULTIPLE BOX A AND B
# A_centers = collect.(Iterators.product(-19.0:0.1:-17.0, 16.0:0.1:18.0))
# A_centers = reshape(A_centers, (length(A_centers),1))
# A_centers = transpose(reduce(hcat, A_centers))

B_centers = collect.(Iterators.product(-98.0:0.1:-92.0, 17.7:0.1:32.0))
B_centers = reshape(B_centers, (length(B_centers),1))
B_centers = transpose(reduce(hcat, B_centers))

corners = [-100, 15, -9, 39]
########################################################################################

function ulamTPT(fname, n_polys, type; h5out = true, extra_suffix = "", rseed = 123)
    # Load data
    data_T = MAT.matopen("$fname.mat")
    x0, xT, y0, yT = read(data_T, "x0", "xT", "y0", "yT")
    x0 = vec(x0)
    y0 = vec(y0)
    xT = vec(xT)
    yT = vec(yT)
    close(data_T)

    # Bin data, get back polygons with polygon centers
    if type == "vor"
        polys, polys_centers = voronoi_binner(x0, y0, xT, yT, n_polys, corners, rseed = rseed) # rseed = 123 is "default"
    elseif type == "hex"
        polys, polys_centers = hexbin_binner(n_polys, corners) 
    elseif type == "reg"        
        polys, polys_centers = square_binner(n_polys, corners) 
    end

    # Ulam's method
    ulam = ulam_nirvana(x0, y0, xT, yT, polys, polys_centers)

    # println(ulam["info"])

    # TPT
    P_closed = ulam["P_closed"]
    P_open = ulam["P_open"]
    pi_open = ulam["pi_open"]
    ALL_inds = 1:length(P_open[:,1])

    ABres = getABinds(ulam["polys_raw"], A_centers, B_centers)
    A_inds = ABres["indsA"]
    B_inds = ABres["indsB"]

    tpt = tpt_infinite_stats(ALL_inds, A_inds, B_inds, P_open, pi_open, P_closed)

    if h5out
        # initialize output file 
        fout_name = "ulamTPT_" * type * "_" * string(n_polys) * extra_suffix * ".h5"
        rm(fout_name, force = true) # manually remove old file if it exists
        fout = h5open(fout_name, "w")

        write_dict_to_h5(fout, "ulam", ulam)
        write_dict_to_h5(fout, "tpt", tpt)

        close(fout)

        println("Output written to: " * fout_name)
    end

    return ulam, tpt
end