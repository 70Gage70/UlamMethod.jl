using NetCDF

include("gdp-process.jl")
include("binner-voronoi.jl")
include("binner-hexbin.jl")
include("binner-square.jl")
include("ulam-nirvana.jl")
include("tpt-infinite.jl")


function tAB_box_T(box_min, box_max, n_box, T_min, T_max, T_step, type; extra_suffix = "")
    boxes = Int64.(floor.(range(box_min, box_max, n_box)))
    deltaT = Int64((T_max - T_min)/T_step) + 1
    Ts = range(T_min, T_max, deltaT)

    n_data_points = length(boxes)*length(Ts)
    data = zeros(n_data_points, 4)

    fout_name = "tAB_box_T" * extra_suffix * ".nc"
    rm(fout_name, force = true) # clear out old file if it exists
    nccreate(fout_name, "tBTdata",  "point" , n_data_points, "data" , 4) # filename, varname, dim1name, dim1, dim2name, dim2 ... 

    north_atlantic_verts = begin [
        -82.5 40 ; -6.06 40 ; -5.9 24.9 ;
        15.3 21.9 ; 15.8 -10 ; -63.7 -10 ; 
        -75.67 5.47 ; -78.24 9.06 ; -80.2 9.65 ; 
        -81.8 9.14 ; -84.2 11.56 ; -84.47 14.12 ; 
        -90.1 15.88 ; -97.67 18.2 ; -100.2 26.08 ; 
        -95.5 31.95 ]
    end

    north_atlantic_edges = begin [
        1 2 ; 2 3 ; 3 4 ; 4 5 ; 5 6 ; 
        6 7 ; 7 8 ; 8 9 ; 9 10 ; 10 11 ;
        11 12 ; 12 13 ; 13 14 ; 14 15 ; 15 16;
        16 1]
    end

    A_centers = [
        -18.0 17.0;
        ]
    
    B_centers = collect.(Iterators.product(-98.0:0.1:-92.0, 17.7:0.1:32.0))
    B_centers = reshape(B_centers, (length(B_centers),1))
    B_centers = transpose(reduce(hcat, B_centers))

    corners = [-100, 15, -9, 39]


    file_name_in = "gdp_raw.csv"
    file_name_in_to_und = "gdp-undrogued.csv"
    file_name_traj = "gdp-all-traj.csv"
    file_name_out_csv = "x0xT_temp.csv"
    file_name_out_mat = "x0xT_temp.mat"

    println("Getting undrogued.")
    rm(file_name_in_to_und, force = true)
    get_undrogued(file_name_in, file_name_in_to_und)

    data_i = 1
    for T in Ts
        println("T = " * string(T))
        println("Generating files for this T.")
        generate_x0xT(file_name_in_to_und, file_name_traj, T)
        select_region(file_name_traj, file_name_out_csv, north_atlantic_verts, north_atlantic_edges)
        csv_to_mat(file_name_out_csv)
        println("Generated.")

        for n_box in boxes
            # Load data
            data_T = MAT.matopen(file_name_out_mat)
            x0, xT, y0, yT = read(data_T, "x0", "xT", "y0", "yT")
            x0 = vec(x0)
            y0 = vec(y0)
            xT = vec(xT)
            yT = vec(yT)
            close(data_T)

            # Bin data, get back polygons with polygon centers
            if type == "vor"
                polys, polys_centers = voronoi_binner(x0, y0, xT, yT, n_box, corners) # rseed = 123 is "default"
            elseif type == "hex"
                polys, polys_centers = hexbin_binner(n_box, corners) 
            elseif type == "reg"        
                polys, polys_centers = square_binner(n_box, corners) 
            end

            # Ulam's method
            ulam = ulam_nirvana(x0, y0, xT, yT, polys, polys_centers, A_centers, B_centers)

            # TPT
            P_open = ulam["P_open"]
            ALL_inds = 1:length(P_open[:,1])
            A_inds = ulam["Ainds"]
            B_inds = ulam["Binds"]
            # P_open = ulam["P_open"]
            pi_open = ulam["pi_open"]
            P_closed = ulam["P_closed"]
            n_box_actual = size(P_open, 1)

            tpt = tpt_infinite_stats(ALL_inds, A_inds, B_inds, P_open, pi_open, P_closed)
    
            data[data_i, :] = [n_box, n_box_actual, T, tpt["tAB"]]
            data_i = data_i + 1
        end
    end

    # clear out temp files
    rm(file_name_in_to_und, force = true)
    rm(file_name_traj, force = true) 
    rm(file_name_out_csv, force = true) 
    rm(file_name_out_mat, force = true)  

    # push to .nc file
    ncwrite(data, fout_name, "tBTdata")
    ncclose(fout_name)

    return "Done."
end

tAB_box_T(20, 1000, 50, 0.5, 5.0, 0.5, "reg", extra_suffix = "_reg_paper")
tAB_box_T(20, 1000, 50, 0.5, 5.0, 0.5, "hex", extra_suffix = "_hex_paper")
tAB_box_T(20, 600, 50, 0.5, 5.0, 0.5, "vor", extra_suffix = "_vor_paper")