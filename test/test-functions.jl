function test_eq(a::Any, b::Any)
    if (length(a) == 0) && (length(b) == 0) && a == b
        return true
    elseif length(a) != length(b)
        return false
    else
        return all(isapprox.(a, b, atol = 1e-13))
    end
end

function ulam_test(ftest::String, ulam_result::UlamResult)
    testf = h5open(ftest)

    # pi_closed
    ulam = ulam_result.pi_closed
    test = read(testf["ulam/pi_closed"])
    if !test_eq(ulam, test)
        @show ulam
        @show test
        @error("pi_closed mismatch")
        return false 
    end 

    # P_closed
    ulam = ulam_result.P_closed
    test = read(testf["ulam/P_closed"])
    if !test_eq(ulam, test)
        @show ulam
        @show test
        @error("P_closed mismatch")
        return false 
    end 

    # polys
    ulam = PolyTable(ulam_result.polys).nodes
    test = read(testf["ulam/polys"])[:,1:2]
    if !test_eq(ulam, test)
        @show ulam
        @show test
        @error("polys mismatch") 
        return false 
    end    

    # polys_dis
    if ulam_result.info.n_polys_dis == 0
        return length(ulam_result.polys_dis) == 0 && length(read(testf["ulam/polys_dis"])) == 0
    else
        ulam = PolyTable(ulam_result.polys_dis).nodes
        test = read(testf["ulam/polys_dis"])[:,1:2]
        if !test_eq(ulam, test)
            @show ulam
            @show test
            @error("polys_dis mismatch")
            return false
        end
    end  

    close(testf)

    return true
end






