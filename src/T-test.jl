using CSV
using Tables
using MAT
using NetCDF

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

function muABplotter(fname, n_polys, type; extra_suffix = "", rseed = 123)
    # Load data
    data_T = MAT.matopen("$fname.mat")
    x0, xT, y0, yT = read(data_T, "x0", "xT", "y0", "yT")
    x0 = vec(x0)
    y0 = vec(y0)
    xT = vec(xT)
    yT = vec(yT)
    close(data_T)

    # fout_name = "TPT_" * type * extra_suffix * ".mat"
    # rm(fout_name, force = true) # need to manually remove old file for some reason
    # fout = matopen(fout_name, "w")
    

    # Bin data, get back polygons with polygon centers
    if type == "vor"
        polys, polys_centers = voronoi_binner(x0, y0, xT, yT, n_polys, corners, rseed = rseed) # rseed = 123 is "default"
        suffix = "vor"
    elseif type == "hex"
        polys, polys_centers = hexbin_binner(n_polys, corners) 
        suffix = "hex"
    elseif type == "reg"        
        polys, polys_centers = square_binner(n_polys, corners) 
        suffix = "reg"
    end

    suffix = suffix * extra_suffix

    # Ulam's method
    ulam = ulam_nirvana(x0, y0, xT, yT, polys, polys_centers, A_centers, B_centers)

    println(ulam["info"])

    # TPT
    P_closed = ulam["P_closed"]
    P_open = ulam["P_open"]
    pi_open = ulam["pi_open"]
    ALL_inds = 1:length(P_open[:,1])
    A_inds = ulam["Ainds"]
    B_inds = ulam["Binds"]
    tpt = tpt_infinite_stats(ALL_inds, A_inds, B_inds, P_open, pi_open, P_closed)

    # Output data
    # write(fout, "P_closed", P_closed)
    # write(fout, "indsA", A_inds)
    # write(fout, "indsB", B_inds)
    # write(fout, "muAB", tpt["muABnorm"])
    # write(fout, "fplus", P_closed)
    # write(fout, "P_closed", tpt["f+"])
    # write(fout, "P_closed", P_closed)
    # write(fout, "P_closed", P_closed)
    # write(fout, "P_closed", P_closed)
    # write(fout, "P_closed", P_closed)
    # write(fout, "P_closed", P_closed)

    csv("indsA_$suffix", A_inds);
    csv("indsB_$suffix", B_inds);
    csv("muAB_$suffix", tpt["muABnorm"]);
    csv("fplus_$suffix", tpt["f+"]);
    csv("trem_$suffix", tpt["t_rem"]);
    csv("polys_$suffix", ulam["polys"]);
    csv("polys_centers_$suffix", ulam["polys_centers"]);
    csv("polys_dis_$suffix", ulam["polys_dis"]);
    csv("pi_open_$suffix", pi_open);
    csv("counts_$suffix", ulam["counts"]);
    csv("leaves_$suffix", ulam["leaves"]);

    println("Output: " * suffix)

    # close(fout)

    return ulam, tpt
end


function tABchecker(fname, n_polys_min, n_polys_max, n_points, type; extra_suffix = "")
    # DATA
    data_T = MAT.matopen("$fname.mat")
    x0, xT, y0, yT = read(data_T, "x0", "xT", "y0", "yT")
    x0 = vec(x0)
    y0 = vec(y0)
    xT = vec(xT)
    yT = vec(yT)
    close(data_T)

    fout_name = "tAB_box_" * string(type) * extra_suffix * ".nc"
    rm(fout_name, force = true) # clear out old file if it exists

    t_res = zeros(n_points, 3)
    n_polys = Int64.(floor.(range(n_polys_min, n_polys_max, n_points)))

    for i = 1:n_points
        println("Percentage complete = " * string(100*i/n_points))
        if type == "vor"
            polys, polys_centers = voronoi_binner(x0, y0, xT, yT, n_polys[i], corners) # rseed = 123
            suffix = "vor"
        elseif type == "hex"
            polys, polys_centers = hexbin_binner(n_polys[i], corners)
            suffix = "hex"
        elseif type == "reg"        
            polys, polys_centers = square_binner(n_polys[i], corners) 
            suffix = "reg"
        end

        # Ulam's method
        ulam = ulam_nirvana(x0, y0, xT, yT, polys, polys_centers, A_centers, B_centers)

        println(ulam["info"])

        # TPT
        P_closed = ulam["P_closed"]
        P_open = ulam["P_open"]
        pi_open = ulam["pi_open"]
        ALL_inds = 1:length(P_open[:,1])
        A_inds = ulam["Ainds"]
        B_inds = ulam["Binds"]
        tpt = tpt_infinite_stats(ALL_inds, A_inds, B_inds, P_open, pi_open, P_closed)

        tAB = tpt["tAB"]
        n_polys_actual = size(P_open, 1)

        t_res[i, :] = [n_polys[i], n_polys_actual, tAB]
    end

    nccreate(fout_name, "data", "trial", n_points, "point", 3)
    ncwrite(t_res, fout_name, "data")
    ncclose(fout_name)

    println("Output is " * fout_name)

    return "Done."

end


function check_centers(fname, ncenters, ntrials; extra_suffix = "")
    # Load data
    data_T = MAT.matopen("$fname.mat")
    x0, xT, y0, yT = read(data_T, "x0", "xT", "y0", "yT")
    x0 = vec(x0)
    y0 = vec(y0)
    xT = vec(xT)
    yT = vec(yT)
    close(data_T)
    
    seeds = zeros(ntrials)
    centers = zeros(ncenters, 2, ntrials)
    tABs = zeros(ntrials)

    fout_name = "CENTERS" * extra_suffix * ".nc"
    rm(fout_name, force = true) # clear out old file if it exists
    nccreate(fout_name, "centers",  "x" , ncenters, "y" , 2, "trial", ntrials) # filename, varname, dim1name, dim2name ... 
    nccreate(fout_name, "seeds",  "seed", ntrials) # filename, varname, dim1name, dim2name ... 
    nccreate(fout_name, "tABs",  "tAB", ntrials)  


    for i = 1:ntrials
        println("")
        println("On trial " * string(i) * " of " * string(ntrials))
        println("")

        rseed = rand(1:ntrials*100)
        seeds[i] = rseed
        polys, polys_centers = voronoi_binner(x0, y0, xT, yT, ncenters, corners, rseed = rseed)

        # Ulam's method
        ulam = ulam_nirvana(x0, y0, xT, yT, polys, polys_centers, A_centers, B_centers)

        println(ulam["info"])

        # TPT
        P_closed = ulam["P_closed"]
        P_open = ulam["P_open"]
        pi_open = ulam["pi_open"]
        ALL_inds = 1:length(P_open[:,1])
        A_inds = ulam["Ainds"]
        B_inds = ulam["Binds"]
        tpt = tpt_infinite_stats(ALL_inds, A_inds, B_inds, P_open, pi_open, P_closed)
        tABs[i] = tpt["tAB"]

        polys_centers = vecvec_to_mat(polys_centers)
        centers[:, :, i] = polys_centers
    end

    ncwrite(centers, fout_name, "centers")
    ncwrite(seeds, fout_name, "seeds")
    ncwrite(tABs, fout_name, "tABs")
    ncclose(fout_name)

    return centers
end


# tABchecker("x0x5-NA-undrogued", 20, 1000, 50, "reg", extra_suffix = "_paper")
# tABchecker("x0x5-NA-undrogued", 20, 1000, 50, "hex", extra_suffix = "_paper")
# tABchecker("x0x5-NA-undrogued", 20, 600, 50, "vor", extra_suffix = "_paper")