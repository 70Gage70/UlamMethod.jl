include(joinpath(@__DIR__, "..", "src/main.jl"))

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

    #JUST MAKE POLYTABLE WORK WITH NO NODES.

    # polys_dis
    ulam = PolyTable(ulam_result.polys_dis).nodes
    test = read(testf["ulam/polys_dis"])[:,1:2]
    if !all(ulam .≈ test) @error("polys_dis mismatch") end  

    close(testf)

    return true
end


function generate_tests(test_cases::Vector)
    for t in test_cases
        traj, domain, fout = t
        res = ulam_method(traj, domain)
        ulam_write(fout, res)
    end

    return
end



