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


function help_smear(xmin, xmax, ymin, ymax; resolution = 0.1)
    centers = collect.(Iterators.product(xmin:resolution:xmax, ymin:resolution:ymax))
    centers = reshape(centers, (length(centers),1))
    centers = transpose(reduce(hcat, centers))
    return centers
end


function ABcorner_standard()
    corners = [-100, 15, -9, 39]

    A_centers = [
        -18.0 17.0;
        ]
    B_centers = help_smear(-98.0, -92.0, 17.7, 32.0, resolution = 0.1)

    return corners, A_centers, B_centers
end

function corners_mid(corners)
    return [(corners[1] + corners[2])/2, (corners[3] + corners[4])/2]
end

# takes in a vector of vectors and converts it to a matrix
function vecvec_to_mat(vvec)
    return reduce(hcat,vvec)'
end

function write_dict_to_h5(fout, group_name, dict)
    g = create_group(fout, group_name)
    for key in collect(keys(dict))
        if dict[key] == []
            # g[key] = HDF5.EmptyArray{Float64}()
            g[key] = zeros(0)
        else
            g[key] = dict[key]
        end
    end
end
