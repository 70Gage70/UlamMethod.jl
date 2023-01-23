"""
Bin the data based a covering of the computational domain by hexagons.

Written by Gage Bonner November 2022
"""
####################################################################################################################################
####################################################################################################################################

"""
Using "odd-r" pointy configuration (top of hexagon is a vertex and every second row is pushed right.) 
Data should be in a rectangle of defined by corners = [xmin, xmax, ymin, ymax].
hex_size is the distance between the center of the hexagon and any of its vertices.
Return dictionary of "nodes" and "edges"
"""

function generate_hexagons(corners, hex_size)
    xmin, xmax, ymin, ymax = corners
    rect_width = xmax - xmin
    rect_height = ymax - ymin
    w = sqrt(3)*hex_size

    n_x = Int64(ceil((rect_width - w/2)/w) + 1)
    n_y = Int64(ceil((rect_height - hex_size)/((3/2)*hex_size)) + 1)

    if (n_x - 1)*w > rect_width # need to trim every second row
        trim_x = true
    else
        trim_x = false
    end

    # extra width and height
    if trim_x
        delta_x = (n_x - 1)*w - rect_width
    else
        delta_x = w/2 + (n_x - 1)*w - rect_width
    end

    delta_y = (n_y - 1)*((3/2)*hex_size) + 2*hex_size - rect_height
    xleft = xmin - delta_x/2
    ytop = ymax - hex_size + delta_y/2

    delta_trim = 0

    if trim_x
        delta_trim = Int64(floor(n_y/2))
    end

    nodes = zeros(6*(n_x*n_y - delta_trim), 2) # x | y
    edges = zeros(Int64, 6*(n_x*n_y - delta_trim), 3) # x | y | cell index
    centers = []
    hex_i = 1
    edge_n = 1
    for i = 1:n_y
        trim_n = 0

        if isodd(i)
            center_x = xleft - w
        else
            center_x = xleft - w + w/2
            if trim_x
                trim_n = 1
            end
        end

        center_y = ytop - (3/2)*hex_size*(i - 1)


        for j = 1:(n_x - trim_n)
            center_x = center_x + w

            push!(centers, [center_x, center_y])
            
            for k = 0:5
                nodes[hex_i + k, :] = [center_x + hex_size*sin((pi/180)*60*k), center_y + hex_size*cos((pi/180)*60*k)]

                if k < 5
                    edges[hex_i + k, :] = [hex_i + k, hex_i + k + 1, edge_n]
                else
                    edges[hex_i + k, :] = [hex_i + k, hex_i, edge_n]
                end
            end

            hex_i = hex_i + 6
            edge_n = edge_n + 1
        end
    end

    return Dict("nodes" => nodes, "edges" => edges, "centers" => centers, "nhex" => edge_n - 1)
end


"""
Clean up for Ulam.
"""

function hexbin_binner(n_polys, corners)
    xmin, xmax, ymin, ymax = corners
    rect_width = xmax - xmin
    rect_height = ymax - ymin
    rect_A = rect_width*rect_height
    hex_size = sqrt(2*rect_A/(3*sqrt(3)*n_polys))

    hexes = generate_hexagons(corners, hex_size)
    polys = []
    nodes = hexes["nodes"]

    for i = 1:6:6*hexes["nhex"]
        this_hex = []
        for k = i:i+5
            push!(this_hex, nodes[k,:])
        end

        push!(this_hex, nodes[i,:]) # close polygon
        push!(polys, this_hex)
    end

    return polys, hexes["centers"]
end