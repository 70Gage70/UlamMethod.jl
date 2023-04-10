function ulam_test(ftest::String, ulam_result::UlamResult)
    testf = h5open(ftest)

    # pi_closed
    ulam = ulam_result.pi_closed
    test = read(testf["ulam/pi_closed"])
    if !all(ulam .≈ test) @error("pi_closed mismatch") end 

    # P_closed
    ulam = ulam_result.P_closed
    test = read(testf["ulam/P_closed"])
    if !all(ulam .≈ test) @error("P_closed mismatch") end 

    # polys
    ulam = PolyTable(ulam_result.polys).nodes
    test = read(testf["ulam/polys"])[:,1:2]
    if !all(ulam .≈ test) @error("polys mismatch") end    

    # polys_dis
    if ulam_result.info.n_polys_dis == 0
        return length(ulam_result.polys_dis) == 0 && length(read(testf["ulam/polys_dis"])) == 0
    else
        ulam = PolyTable(ulam_result.polys_dis).nodes
        test = read(testf["ulam/polys_dis"])[:,1:2]
        if !all(ulam .≈ test) @error("polys_dis mismatch") end
    end  

    close(testf)

    return true
end






