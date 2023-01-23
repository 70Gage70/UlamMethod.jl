"""
Bin the data based a covering of the computational domain by squares.

Written by Gage Bonner November 2022
"""
####################################################################################################################################
####################################################################################################################################

function square_binner(n_polys, corners)
    xmin, xmax, ymin, ymax = corners
    w = xmax - xmin
    L = ymax - ymin
    s = sqrt(w*L/n_polys)

    n_x = Int64(ceil(w/s))
    n_y = Int64(ceil(L/s))

    delta_x = n_x*s - w
    delta_y = n_y*s - L

    y_top = ymax + delta_y/2
    x_left = xmin - delta_x/2

    polys = []
    centers = []

    for i = 1:n_y
        x_left = xmin - delta_x/2

        for j = 1:n_x
            poly = [[x_left, y_top - s], [x_left, y_top], [x_left + s, y_top], [x_left + s, y_top - s], [x_left, y_top - s]] # note polygon is closed
            push!(polys, poly)

            x_center = x_left + s/2
            y_center = y_top - s/2
            push!(centers, [x_center, y_center])

            x_left = x_left + s
        end

        y_top = y_top - s
    end

    return polys, centers
end
