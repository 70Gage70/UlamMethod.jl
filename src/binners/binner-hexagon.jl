"""
    binner_hexagon(domain)

Cover the rectangle defined by `domain.corners` by a uniform grid of hexagons.

Hexagons near the edges may overlap with the rectangle.
"""
function binner_hexagon(domain::UlamDomain)
    xmin, xmax, ymin, ymax = domain.corners
    n_polys = domain.poly_number

    # The user requests a number of hexagons, we have to then calculate the approximate size by ratios of areas
    rect_width = xmax - xmin
    rect_height = ymax - ymin
    rect_A = rect_width*rect_height
    hex_size = sqrt(2*rect_A/(3*sqrt(3)*n_polys))

    # w is the width of a hexagon, n_x and n_y are the numbers of hexagons in each direction
    w = sqrt(3)*hex_size
    n_x = Int64(ceil((rect_width - w/2)/w) + 1)
    n_y = Int64(ceil((rect_height - hex_size)/((3/2)*hex_size)) + 1)

    # trim_x takes into account whether there will be one too many hexagons outside the given 
    # rectangle, which happens exactly when (n_x - 1)*w > rect_width
    # in this case, we trim every second row (we don't have to trim every row since they are offset)
    # and this has to be taken into account for the entire calculation
    if (n_x - 1)*w > rect_width
        trim_x = true
    else
        trim_x = false
    end

    # delta_x and delta_y ensure the covering is tight
    if trim_x
        delta_x = (n_x - 1)*w - rect_width
    else
        delta_x = w/2 + (n_x - 1)*w - rect_width
    end

    delta_y = (n_y - 1)*((3/2)*hex_size) + 2*hex_size - rect_height

    # we build up the covering in rows from left to right starting from the top
    xleft = xmin - delta_x/2
    ytop = ymax - hex_size + delta_y/2

    delta_trim = 0
    if trim_x
        delta_trim = Int64(floor(n_y/2))
    end

    nodes = zeros(6*(n_x*n_y - delta_trim), 2) # x | y
    hex_i = 1

    for i = 1:n_y
        # first we calculate the center of the leftmost hexagon we're creating

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

        # now we look across the rows; note the one less hexagon when a trim is required
        for j = 1:(n_x - trim_n)
            center_x = center_x + w
            
            # we construct the vertices of the current hexagon
            # vertices of a hexagon are 60 degrees apart with respect to their center
            for k = 0:5
                nodes[hex_i + k, :] = [center_x + hex_size*sin((pi/180)*60*k), center_y + hex_size*cos((pi/180)*60*k)]
            end

            hex_i = hex_i + 6
        end
    end

    polys = [UlamPolygon(nodes[i:i+5,:]) for i = 1:6:size(nodes, 1)]

    return polys
end

# provides an interface for binner-triangle
function binner_hexagon(xmin::Real, xmax::Real, ymin::Real, ymax::Real, poly_number::Integer)
    domain = UlamDomain(xmin, xmax, ymin, ymax, poly_number = poly_number)

    return binner_hexagon(domain)
end