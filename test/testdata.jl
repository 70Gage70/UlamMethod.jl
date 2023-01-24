include(joinpath(@__DIR__, "..", "main.jl"))

function ulam_method_test(f_in, n_polys, type, corners, f_standard)
    ulam = ulam_method(f_in, n_polys, type, corners, h5out = false)
    fid = h5open(f_standard, "r")

    for key in collect(keys(fid["ulam"]))
        if !(key in collect(keys(ulam)))
            return "No ulam key: " * string(key)
        else
            expected = read(fid["ulam"][key])
            calculated = ulam[key]
            if (typeof(expected) == String) && (typeof(calculated) == String)
                equal = expected == calculated
            else
                equal = all(expected .≈ calculated)
            end

            if !equal
                return "Bad ulam key: " * string(key) * ". Expected " * string(expected) * " but calculated " * string(calculated)
            end
        end
    end

    close(fid)

    return true
end


function ulam_method_test(ulam, A_centers, B_centers, f_standard)
    tpt = tpt_from_ulam(ulam, A_centers, B_centers, h5out = false)
    fid = h5open(f_standard, "r")

    for key in collect(keys(fid["tpt"]))
        if !(key in collect(keys(tpt)))
            return "No tpt key: " * string(key)
        else
            expected = read(fid["tpt"][key])
            calculated = tpt[key]
            if (typeof(expected) == String) && (typeof(calculated) == String)
                equal = expected == calculated
            else
                equal = all(expected .≈ calculated)
            end

            if !equal
                return "Bad tpt key: " * string(key) * ". Expected " * string(expected) * " but calculated " * string(calculated)
            end
        end
    end

    close(fid)

    return true
end

# Computational domain
corners = [-100, 15, -9, 39]

# AB Locations
A_centers = [
    -18.0 17.0;
    ]
B_centers = AB_smear(-98.0, -92.0, 17.7, 32.0, resolution = 0.1)
