using MAT
using HDF5

function process_infile(infile::String)
    extension = infile[findlast(==('.'), infile)+1:end]
    if extension == "mat"
        data_T = MAT.matopen(infile)
        x0, xT, y0, yT = read(data_T, "x0", "xT", "y0", "yT")
        close(data_T)
    elseif extension == "h5"
        data_T = HDF5.h5open(infile)
        x0, xT, y0, yT = read(data_T["x0"]), read(data_T["xT"]), read(data_T["y0"]), read(data_T["yT"])
        close(data_T)
    else
        error("File type not supported.") 
    end

    x0 = vec(x0)
    y0 = vec(y0)
    xT = vec(xT)
    yT = vec(yT)

    return x0, y0, xT, yT
end