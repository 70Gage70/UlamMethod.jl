using MAT 
using DelimitedFiles

include("earth-polygons.jl")
include("tpt-infinite.jl")
# include("ulam-regular.jl")
include("ulam-regular-opt.jl")

###########################
###### Load in Data
###########################

# gdp_data = matopen("gdp_final.mat")
gdp_data = matopen("x0x5-NA-undrogued.mat")
x0, xT, y0, yT = read(gdp_data, "x0", "xT", "y0", "yT")
close(gdp_data)

###########################
###### Ulam's Method
###########################

A_center = [-18, 17]
B_center = [-95, 27]
corners = [-100, 15, -9, 39]

divs = exp10.(range(log10(0.6), stop=log10(10), length=30))
div_box = zeros(length(divs), 2)

for i = 1:length(divs)
    div = divs[i]
    display("Percentage complete = " * string(i/length(divs)))
    @time ulamRES = ulamreg(vec(x0), vec(y0), vec(xT), vec(yT), div, A_center, B_center, corners)
    display(ulamRES["info"])

    ###########################
    ###### TPT
    ###########################

    P_open = ulamRES["P_open"] 
    pi_open = ulamRES["pi_open"] 
    
    ALL_inds = 1:length(P_open[:, 1])
    A_inds = ulamRES["Ainds"] 
    B_inds = ulamRES["Binds"]
    gps = ulamRES["gps"]
    
    tpt = tpt_infinite_stats(ALL_inds, A_inds, B_inds, P_open, pi_open) 

    div_box[i, :] = [ulamRES["div_actual"], tpt["tAB"]*5/365]
end

writedlm( "div_box2.csv",  div_box, ',')
