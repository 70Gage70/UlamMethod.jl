using HDF5

include("binner-voronoi.jl")
include("binner-hexbin.jl")
include("binner-square.jl")

include("ulam-nirvana.jl")
include("tpt-infinite.jl")

include("helpers.jl")
include("earth-polygons.jl")
include("process_infile.jl")

#########################################################################################################

function ulam_method(fname, n_polys, bin_type, corners; sto_type = "data", source_centers = corners_mid(corners), h5out = true, extra_suffix = "", rseed = 123)
    x0, y0, xT, yT = process_infile(fname)

    if bin_type == "vor"
        polys, polys_centers = voronoi_binner(x0, y0, xT, yT, n_polys, corners, rseed = rseed) # rseed = 123 is "default"
    elseif bin_type == "hex"
        polys, polys_centers = hexbin_binner(n_polys, corners) 
    elseif bin_type == "reg"        
        polys, polys_centers = square_binner(n_polys, corners) 
    end

    ulam = ulam_nirvana(x0, y0, xT, yT, polys, polys_centers, sto_type = sto_type, source_centers = source_centers)

    if h5out
        # initialize output file 
        fout_name = "ulam_" * bin_type * "_" * string(n_polys) * extra_suffix * ".h5"
        rm(fout_name, force = true) # manually remove old file if it exists
        fout = h5open(fout_name, "w")
        write_dict_to_h5(fout, "ulam", ulam)
        close(fout)

        println("Output written to: " * fout_name)
    end

    return ulam
end


function AB_smear(xmin, xmax, ymin, ymax; resolution = 0.1)
    centers = collect.(Iterators.product(xmin:resolution:xmax, ymin:resolution:ymax))
    centers = reshape(centers, (length(centers),1))
    centers = transpose(reduce(hcat, centers))
    return centers
end


function tpt_from_ulam(ulam, A_centers, B_centers; h5out = true, extra_suffix = "")
    P_closed = ulam["P_closed"]
    P_open = ulam["P_open"]
    pi_open = ulam["pi_open"]
    ALL_inds = 1:length(P_open[:,1])

    A_inds = getABinds(ulam["polys"], A_centers)
    B_inds = getABinds(ulam["polys"], B_centers)

    tpt = tpt_infinite_stats(ALL_inds, A_inds, B_inds, P_open, pi_open, P_closed)

    if h5out
        # initialize output file 
        fout_name = "TPT" * extra_suffix * ".h5"
        rm(fout_name, force = true) # manually remove old file if it exists
        fout = h5open(fout_name, "w")
        write_dict_to_h5(fout, "tpt", tpt)
        close(fout)

        println("Output written to: " * fout_name)
    end

    return tpt
end