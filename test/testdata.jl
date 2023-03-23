include(joinpath(@__DIR__, "..", "src/main.jl"))

function ulam_method_tpt_test(f_in, n_polys, type, corners, A_centers, B_centers, f_ulam, f_tpt)
    ulam = ulam_method(f_in, n_polys, type, corners, h5out = false)
    fid = h5open(f_ulam, "r")

    for key in collect(keys(ulam))
        if !(key in collect(keys(fid["ulam"])))
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

    tpt = tpt_from_ulam(ulam, A_centers, B_centers, h5out = false)
    fid = h5open(f_tpt, "r")

    for key in collect(keys(tpt))
        if !(key in collect(keys(fid["tpt"])))
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

# GENERATING TEST REFERENCES

# Computational domain
corners = [-100, 15, -9, 39];

# AB Locations
A_centers = [
    -18.0 17.0;
    ];
B_centers = AB_smear(-98.0, -92.0, 17.7, 32.0, resolution = 0.1);


test_cases = [["reg", 500], ["hex", 500], ["vor", 50]]
test_data_file = "x0x5-NA-undrogued.mat"

function generate_tests(test_data_file, test_cases)
    for i = 1:length(test_cases)
        type = test_cases[i][1]
        n_polys = test_cases[i][2]
        ulam = ulam_method(test_data_file, n_polys, type, corners)
        tpt = tpt_from_ulam(ulam, A_centers, B_centers, extra_suffix = "_" * string(type) * "_" * string(n_polys))
    end

    return "Tests generated."
end

