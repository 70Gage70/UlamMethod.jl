"""
General helper functions.
"""
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

import MAT
import HDF5

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

# takes in a vector of vectors and converts it to a matrix
function vecvec_to_mat(vvec)
    return reduce(hcat,vvec)'
end

function write_dict_to_h5(fout, group_name, dict)
    g = create_group(fout, group_name)
    # println(collect(keys(dict)))
    for key in collect(keys(dict))
        # println(key)
        # println(dict[key])
        if key != "polys_raw"
            if dict[key] == []
                # g[key] = HDF5.EmptyArray{Float64}()
                g[key] = zeros(0)
            else
                g[key] = dict[key]
            end
        end
    end
end
