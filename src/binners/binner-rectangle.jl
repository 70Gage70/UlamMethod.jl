"""
    binner_rectangle(domain)

Cover the rectangle defined by `domain.corners` by a tight uniform grid of rectangles.

The rectangles are chosen to be as close to squares as possible with the constraint that they must exactly fit `domain.corners`.
"""
function binner_rectangle(domain::UlamDomain)
    n_polys = domain.poly_number
    xmin, xmax, ymin, ymax = domain.corners

    w = xmax - xmin 
    L = ymax - ymin
    n_x = sqrt(n_polys/(L/w))
    n_y = (L/w)*n_x

    n_x = round(Int64, n_x)
    n_y = round(Int64, n_y)

    # construct a grid of vertices
    points = Iterators.product(
        range(xmin, xmax, length = n_x + 1), # + 1 since N points is N - 1 squares
        range(ymin, ymax, length = n_y + 1)
    )

    points = collect.(points) # convert from tuples to matrix of vectors

    polys = Vector{UlamPolygon{Float64}}()
    for i = 1:size(points, 1) - 1
        for j = 1:size(points, 2) - 1
            points4 = reshape(points[i:i + 1, j:j + 1], 4, 1)[[1, 3, 4, 2]] # a vector of vectors of ordered vertices of the rectange
            verts = Matrix(hcat(points4...)') # converting to a 4 x 2 matrix of points

            push!(polys, UlamPolygon(verts))
        end
    end

    return polys
end
