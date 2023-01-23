"""
General veriables and data for Ulam/TPT computations.
"""
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

# import CSV
# import Tables
import MAT
import HDF5

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

# DRIFTER DATA
# gdp_data = MAT.matopen("x0x5-NA-undrogued.mat")
# x0, xT, y0, yT = read(gdp_data, "x0", "xT", "y0", "yT")
# x0 = vec(x0)
# y0 = vec(y0)
# xT = vec(xT)
# yT = vec(yT)
# close(gdp_data)

# SINGLE BOX A AND B
# A_centers = [
#     -18.0 17.0;
#     ]
# B_centers = [
#     -95.0 27.0;
#     ]

# MULTIPLE BOX A AND B
# A_centers = collect.(Iterators.product(-19.0:0.1:-17.0, 16.0:0.1:18.0))
# A_centers = reshape(A_centers, (length(A_centers),1))
# A_centers = transpose(reduce(hcat, A_centers))

# B_centers = collect.(Iterators.product(-98.0:0.1:-92.0, 17.7:0.1:32.0))
# B_centers = reshape(B_centers, (length(B_centers),1))
# B_centers = transpose(reduce(hcat, B_centers))

# corners = [-100, 15, -9, 39]

# USEFUL FUNCTIONS
# function csv(name, arr)
#     CSV.write("$name.csv", Tables.table(arr), delim = ',', writeheader = false)
#     return
# end

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
