"""
Data is downloaded from http://osmc.noaa.gov/erddap/tabledap/gdp_interpolated_drifter.html with the following elements selected. All others are unselected.

- WM0, at least 1000000 (integer)
- longitude (degrees_east)
- latitude (degrees_north)
- time (UTC)
- drogue_lost_date (UTC) [set the constraints to be between min and max time manually, 
    otherwise this will include data whose drogue is never lost, or for which the status of the drogue is unknown]

Sort data by WMO, then by time.

File is saved as gdp_raw.csv

Written by Gage Bonner September 2022
"""

###################################################################################################
###################################################################################################
###################################################################################################

using Dates # For handling UTC time
using PolygonInbounds # used for algorithm which determines if a point is inside a polygon
using MAT # used for outputting a mat file if needed


"""
Takes in a line from readline(file.csv) and returns a dict of the values.
"""

function parse_gdp_line(readline)

    line = split(readline, ",")

    id = parse(Int64, line[1]) # drifter ID
    lon = parse(Float64, line[2]) # drifter longitude (between -180 and 180)
    lat = parse(Float64, line[3]) # drifter latitude (between -90 and 90)               
    time = DateTime(line[4][1:end - 1]) # UTC time, remove last character because DateTime doesn't want the trailing "Z"
    d_lost = DateTime(line[5][1:end - 1]) # UTC time drifter is lost

    returndict = Dict(
    "id" => id,
    "lon" => lon,
    "lat" => lat,
    "time" => time,
    "d_lost" => d_lost)

    return returndict    
end

"""
Takes raw 6h interpolated GDP file and returns a CSV file with four columns
|ID|time|lon|lat
Time measured in days since drogue was lost.
Time increments of trajectories are all 0.25 days (6 hours); if an increment is longer, 
a new trajectory is started.
"""

function get_undrogued(file_name_in, file_name_out)
    raw_data = open(file_name_in)
    cleaned_data = open(file_name_out,"w")
    
    readline(raw_data) # variables
    readline(raw_data) # units
    
    current = parse_gdp_line(readline(raw_data))
    
    last_id = current["id"]
    d_lost = current["d_lost"]
    last_t = current["time"]
    n_traj = 1
    double_obs = 1 
    while !eof(raw_data) # while we haven't reached the end of the raw_data file
        next = parse_gdp_line(readline(raw_data))
        
        
        if (next["id"] != last_id) || (next["d_lost"] != d_lost) # we're on to a new trajectory 
            d_lost = next["d_lost"]
            last_id = next["id"]
            n_traj = n_traj + 1
        else # we're still on the same trajectory
            if next["time"] >= d_lost # we're on or at the time the drogue is lost

                dt = Dates.Hour(next["time"] - last_t)/Dates.Hour(24) # days since last observation on this trajectory            
                
                if dt > 0.25 # start a new trajectory, we're too far away (temporally) from the last observation
                    d_lost = next["d_lost"] # shouldn't actually change
                    last_id = next["id"] # shouldn't actually change
                    n_traj = n_traj + 1               
                elseif dt == 0.25
                    traj_id = string(n_traj)
                    delta_t = string(Dates.Hour(next["time"] - d_lost)/Dates.Hour(24)) # Days since drogue lost
                    lon = string(next["lon"])
                    lat = string(next["lat"])
                    write(cleaned_data, traj_id * "," * delta_t * "," * lon * "," * lat * "\n")
                else    
                    double_obs = double_obs + 1 # for some reason, sometimes there are two observations at the same time; ignore
                end 
            end
        end
        
        last_t = next["time"] # time of most recent observation
    end

    close(cleaned_data)
    close(raw_data)
    
#     display(double_obs)
    return "Completed."
end 

"""
Takes in a line from the output of get_undrogued and returns a dict of the values.
"""

function parse_cleaned_gdp_line(readline)

    line = split(readline, ",")

    id = parse(Int64, line[1]) # trajectory id
    time = parse(Float64, line[2]) # days since lost drogue
    lon = parse(Float64, line[3]) # longitude             
    lat = parse(Float64, line[4]) # latitude

    returndict = Dict(
    "id" => id,
    "time" => time,
    "lon" => lon,
    "lat" => lat)

    return returndict    
end

"""
Takes in a cleaned 6h GDP file from get_undrogued and generates initial and final observation points
Input T is the the time in days between initial and final obs, generally should be 2 - 5 days.
"""

function generate_x0xT(file_name_in, file_name_out, T)
    T = Int64(4*T) # because e.g. 5 days is 4*5 = 20 increments of 6 hours
    raw_data = open(file_name_in)
    cleaned_data = open(file_name_out,"w")

    current = parse_cleaned_gdp_line(readline(raw_data))
    last_id = current["id"]
    this_traj_lon = Vector{Float64}()
    this_traj_lat = Vector{Float64}()   
    
    while !eof(raw_data) 
        next = parse_cleaned_gdp_line(readline(raw_data))
        
        if next["id"] == last_id # we're on the same trajectory
            push!(this_traj_lon, next["lon"])
            push!(this_traj_lat, next["lat"])
        else # switched trajectories, collect and write data if possible
            if length(this_traj_lon) >= T # otherwise trajectory is too short
                for i = 1:length(this_traj_lon) - (T - 1)
                    x0 = string(this_traj_lon[i])
                    y0 = string(this_traj_lat[i])
                    xT = string(this_traj_lon[i + T - 1])
                    yT = string(this_traj_lat[i + T - 1])
                    write(cleaned_data, x0 * "," * y0 * "," * xT * "," * yT * "\n")
                end
            end
            
            last_id = next["id"]
            this_traj_lon = Vector{Float64}()
            this_traj_lat = Vector{Float64}() 
            
        end
    end
                
    close(cleaned_data)
    close(raw_data)      
        
    return "Completed"
        
end

"""
picks out points from generate_x0xT trajectories which lie in a given polygon
"""

function select_region(file_name_in, file_name_out, verts, edges)
    raw_data = open(file_name_in)
    cleaned_data = open(file_name_out,"w")
    
    while !eof(raw_data)
        line = split(readline(raw_data), ",")

        x0 = parse(Float64, line[1]) 
        y0 = parse(Float64, line[2]) 
        xT = parse(Float64, line[3])           
        yT = parse(Float64, line[4]) 

        p1 = inpoly2([x0, y0], verts, edges)[1,1]
        p2 = inpoly2([xT, yT], verts, edges)[1,1]

        if p1 && p2
            write(cleaned_data, string(x0) * "," * string(y0) * "," * string(xT) * "," * string(yT) * "\n")
        end
    end

    close(cleaned_data)
    close(raw_data)     
    
    return "Completed"
end

"""
Converts the csv output from select_region to a mat file
"""

function csv_to_mat(file_name_in)
    raw_data = open(file_name_in)
    data = readlines(raw_data)
    x0 = Vector{Float64}()
    y0 = Vector{Float64}()
    xT = Vector{Float64}()
    yT = Vector{Float64}()
    
    for point in data
        point = split(point, ",")
        push!(x0, parse(Float64, point[1]))
        push!(y0, parse(Float64, point[2]))
        push!(xT, parse(Float64, point[3]))
        push!(yT, parse(Float64, point[4]))
    end
        
    file = matopen(file_name_in[1:end-3] * "mat", "w")
    write(file, "x0", x0)
    write(file, "y0", y0)
    write(file, "xT", xT)
    write(file, "yT", yT)
    
    close(file)
    close(raw_data)
end 

"""
Inputs:
file_name_in: name of the GDP file (explained in header)
file_name_out: name of final output (something like ____.csv probably)
T: time in days between initial and final observations
verts/coords: vertices and coordinates of polycon where data is to be constrained in

Outputs:
gdp-undrogued.csv, csv of observations with only trajectories which start from the undrogued started
gdp-all-traj.csv, csv of all x0xT pairs for time T
file_name_out.csv, csv of final trajectories constrained to be in the supplied region
file_name_out.mat, matlab file of final trajectories constrained to be in the supplied region
"""

function gdp_to_ulam_undrogued(file_name_in, file_name_out, T, verts, edges)
    get_undrogued(file_name_in, "gdp-undrogued.csv")
    display("Removed drogued trajectories.")
    generate_x0xT("gdp-undrogued.csv", "gdp-all-traj.csv", T)
    display("Generated trajectories over the entire dataset.")
    select_region("gdp-all-traj.csv", file_name_out, verts, edges)
    display("Generated trajectories restricted to the given area.")
    csv_to_mat(file_name_out)
    display("Generated csv and matlab file.")
    display("Done.")
    return true
end