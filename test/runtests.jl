using UlamMethod
using Test
using HDF5

function ulamTPTtest(fin, n_polys, type, f_standard)
    ulam, tpt = ulamTPT(fin, n_polys, type, h5out = false)
    fid = h5open(f_standard, "r")

    # test Ulam's method
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

    # test TPT
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

inpath = joinpath(@__DIR__, "x0x5-NA-undrogued")
regpath = joinpath(@__DIR__, "ulamTPT_reg_500.h5")
hexpath = joinpath(@__DIR__, "ulamTPT_hex_500.h5")
vorpath = joinpath(@__DIR__, "ulamTPT_vor_50.h5")

@testset "UlamMethod.jl" begin
    @test ulamTPTtest(inpath, 500, "reg", regpath) == true
    @test ulamTPTtest(inpath, 500, "hex", hexpath) == true
    @test ulamTPTtest(inpath, 50, "vor", vorpath) == true
end