using MAT 
using DelimitedFiles

include("earth-polygons.jl")
include("tpt-infinite.jl")
include("ulam-bisection.jl")

###########################
###### Load in Data
###########################

# gdp_data = matopen("gdp_final.mat") # Gage's processing
gdp_data = matopen("x0x5-NA-undrogued.mat")
x0, xT, y0, yT = read(gdp_data, "x0", "xT", "y0", "yT")
x0 = vec(x0)
y0 = vec(y0)
xT = vec(xT)
yT = vec(yT)
close(gdp_data)

###########################
###### Ulam's Method
###########################

n_boxes = 100
A_center = [-18, 17]
B_center = [-95, 27]
corners = [-100, 15, -9, 39]

# ###########################
# ###### TPT Infinite
# ###########################

n_boxes = exp10.(range(log10(100), stop=log10(10000), length=30))
div_box = zeros(length(divs), 2)

for i = 1:length(divs)
    div = divs[i]
    display("Percentage complete = " * string(i/length(divs)))
    @time ulamRESbi = ulam_bisection(x0, xT, y0, yT, corners, n_boxes[i], A_center, B_center)
    display(ulamRESbi["info"])

    ###########################
    ###### TPT
    ###########################

    P_open = ulamRESbi["P_open"] 
    pi_open = ulamRESbi["pi_open"] 
    
    ALL_inds = 1:length(P_open[:, 1])
    A_inds = ulamRESbi["Ainds"] 
    B_inds = ulamRESbi["Binds"]
    gps = ulamRESbi["gridboxes"]
    
    tpt = tpt_infinite_stats(ALL_inds, A_inds, B_inds, P_open, pi_open) 

    div_box[i, :] = [n_boxes[i], tpt["tAB"]*5/365]
end

writedlm( "div_box_bisection.csv",  div_box, ',')
