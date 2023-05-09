"""
    binner_triangle(domain)

Cover the rectangle defined by `domain.corners` by a uniform grid of equilateral triangles.

Triangles near the edges may overlap with the rectangle.
"""
function binner_triangle(domain::UlamDomain)
    xmin, xmax, ymin, ymax = domain.corners
    poly_number = ceil(Int64, domain.poly_number/6) # each hexagon is the union of 6 triangles

    hexagons = binner_hexagon(xmin, xmax, ymin, ymax, poly_number)
    polys = Vector{UlamPolygon{Float64}}()

    for hex in hexagons
        for tri_inds in [[1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 1]]
            tri_nodes = vcat(hex.nodes[tri_inds,:], hex.center)
            push!(polys, UlamPolygon(tri_nodes))
        end
    end

    return polys
end