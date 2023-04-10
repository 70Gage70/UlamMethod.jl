ftest = joinpath(@__DIR__, "x0x5-NA-undrogued.h5")

t1 = [
    UlamTrajectories(ftest), 
    UlamDomain([-100, 15, -9, 39]..., poly_type = "sqr", poly_number = 500),
    joinpath(@__DIR__, "ulam_sqr_500.h5")
    
]

t2 = [
    UlamTrajectories(ftest), 
    UlamDomain([-100, 15, -9, 39]..., poly_type = "hex", poly_number = 500),
    joinpath(@__DIR__, "ulam_hex_500.h5")]

t3 = [
    UlamTrajectories(ftest), 
    UlamDomain([-100, 15, -9, 39]..., poly_type = "vor", poly_number = 50),
    joinpath(@__DIR__, "ulam_vor_50.h5")
]

test_cases = [
    t1,
    t2,
    t3
]