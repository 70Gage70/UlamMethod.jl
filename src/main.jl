using HDF5
# using .UlamTypes

include("helpers.jl")

include("ulam-nirvana.jl")
include("ulam-binner.jl")

include("earth-polygons.jl")

#########################################################################################################


"""
    ulam_method(traj, domain)

Run the high-level Ulam method and return an [`UlamResult`](@ref).

### Arguments
- `traj`: An [`UlamTrajectories`](@ref); contains the trajectory data.
- `domain`: An [`UlamDomain`](@ref); contains the domain specification.
"""
function ulam_method(traj::UlamTrajectories, domain::UlamDomain)
    polys = ulam_binner(traj, domain)
    ulam = ulam_nirvana(traj, domain, polys)

    return ulam
end

"""
    ulam_write(outfile, ulam_result; dir_name, overwrite, P_out)

Write `ulam_result` to the file `outfile`, which must be in the `.h5` format.

If `outfile` does not exist, it will be created. Results are written to the directory specified by `dir_name`.

### Optional Arguments
- `dir_name`: The name of the directory that `ulam_result` is written to, default `"ulam"`. Directories can be nested, e.g. `"trial1/ulam"`.
- `overwrite`: If `true`, the directory `dir_name` will overwrite a directory with the same name if it exists in `outfile`. Default `false`.
- `P_out`: If `false`, `P_closed` is not written to file. Default `true`.
"""
function ulam_write(
    outfile::String, 
    ulam_result::UlamResult; 
    dir_name::String = "ulam", 
    overwrite::Bool = false,
    P_out::Bool = true)

    @assert outfile[end-2:end] == ".h5" "The output file must be of the form filename.h5"
    fout = h5open(outfile, "cw")

    if dir_name in keys(fout)
        if !overwrite
            @assert !(dir_name in keys(fout)) "This file already has a group with the name: $(dir_name). Pass `overwrite = true` to force a replacement."
        end

        delete_object(fout, dir_name)
    end

    g = create_group(fout, dir_name)

    g["n_polys"] = length(ulam_result.polys)
    g["n_polys_dis"] = length(ulam_result.polys_dis)
    g["pi_closed"] = ulam_result.pi_closed
    g["counts"] = ulam_result.counts

    # The polygons are output in to an n_polys x 3 matrix. The first two columns of the matrix are
    # the (x, y) nodes and the third column is the index of the polygon that node belongs to.
    g["polys"] = [PolyTable(ulam_result.polys).nodes ;; PolyTable(ulam_result.polys).edges[:,3]]

    # Similarly for disconnected polygons, but we explicitly write an empty array if there are no
    # disconnected polygons.
    g["polys_dis"] = ulam_result.info.n_polys_dis == 0 ? zeros(0) : 
    [PolyTable(ulam_result.polys_dis).nodes ;; PolyTable(ulam_result.polys_dis).edges[:,3]]

    if P_out g["P_closed"] = ulam_result.P_closed end

    close(fout)

    @info "UlamResult written to $(outfile)."

    return
end