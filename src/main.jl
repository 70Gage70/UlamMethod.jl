using HDF5
using .UlamTypes

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
    ulam_write(outfile, ulam_result; P_out)

Write `ulam_result` to the file `outfile`, which must be in the `.h5` format. 

### Optional Arguments
- `P_out`: If false, `P_closed` is not written to file since it can sometimes be very large. Default true.
"""
function ulam_write(outfile::String, ulam_result::UlamResult; P_out::Bool = true)
    @assert outfile[end-2:end] == ".h5" "The output file must be of the form filename.h5"

    fout = h5open(outfile, "w")
    g = create_group(fout, "ulam")

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