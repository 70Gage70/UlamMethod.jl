"""
    square_binner(ulam_problem)

Cover the computational domain in `ulam_problem` by a uniform grid of squares and return an [`UlamCovering`](@ref).
"""
function square_binner(traj::UlamTrajectories, domain::UlamDomain)::UlamCovering
    n_polys = domain.bin_number
    xmin, xmax, ymin, ymax = domain.corners
    w = xmax - xmin
    L = ymax - ymin
    s = sqrt(w*L/n_polys)

    n_x = Int64(ceil(w/s))
    n_y = Int64(ceil(L/s))

    delta_x = n_x*s - w
    delta_y = n_y*s - L

    y_top = ymax + delta_y/2
    x_left = xmin - delta_x/2

    polys = Vector{UlamPolygon}()

    for i = 1:n_y
        x_left = xmin - delta_x/2

        for j = 1:n_x
            poly = [
                x_left y_top - s;
                x_left y_top;
                x_left + s y_top;
                x_left + s y_top - s
                ]

            push!(polys, UlamPolygon(poly))

            x_left = x_left + s
        end

        y_top = y_top - s
    end

    return UlamCovering(polys)
end
