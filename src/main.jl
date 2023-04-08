using HDF5

include("helpers.jl")

include("ulam-nirvana.jl")
include("ulam-binner.jl")

include("earth-polygons.jl")

#########################################################################################################


"""
    ulam_method(traj, domain)

The high-level Ulam method; cover the domain by polygons and then construct the transition probability matrix.
"""
function ulam_method(traj::UlamTrajectories, domain::UlamDomain)
    polys = ulam_binner(traj, domain)
    ulam = ulam_nirvana(traj, domain, polys)

    return ulam
end

# function ulam_method(fulam::String)::UlamResult
#     traj = UlamTrajectories(fulam)
#     domain = UlamDomain(fulam)
#     return ulam_method(traj, domain)
# end

# need to write a function which writes an UlamResult to a file

# if h5out
#     # initialize output file 
#     fout_name = "ulam_" * bin_type * "_" * string(n_polys) * extra_suffix * ".h5"
#     rm(fout_name, force = true) # manually remove old file if it exists
#     fout = h5open(fout_name, "w")
#     write_dict_to_h5(fout, "ulam", ulam)
#     close(fout)

#     println("Output written to: " * fout_name)
# end





