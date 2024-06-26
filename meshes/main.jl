using Meshes, CairoMakie
# coords(point).x to get values

# 1d boundary
boundary1d = Segment((0,), (1,))

# 2d boundary
boundary2d = Ngon([(0,0),(6,0),(1,7),(1,6)]...)

abstract type AbstractBinner{Dim} end

struct LineBinner <: AbstractBinner{1}
    nbins::Int64
end

struct RectangleBinner <: AbstractBinner{2}
    nbins::Int64
end

function bin(boundary::Chain, binner::LineBinner)
    nbins = binner.nbins
    bbox = Meshes.boundingbox(boundary)
    grid = CartesianGrid(bbox.min, bbox.max, dims = (nbins, ))
    return GeometrySet([elem for elem in elements(grid) if intersects(elem, boundary)])
end 

function bin(boundary::Polygon, binner::RectangleBinner)
    nbins = binner.nbins

    bbox = Meshes.boundingbox(boundary)
    bbox_W = coords(bbox.max).x - coords(bbox.min).x
    bbox_L = coords(bbox.max).y - coords(bbox.min).y
    
    n_x = sqrt(bbox_W*nbins/bbox_L) |> x -> round(Int64, x)
    n_y = sqrt(bbox_L*nbins/bbox_W) |> x -> round(Int64, x)
    
    grid = CartesianGrid(bbox.min, bbox.max, dims = (n_x, n_y))
    return GeometrySet([elem for elem in elements(grid) if intersects(elem, boundary)])
end 


test_bin1d = bin(boundary1d, LineBinner(500))
test_bin2d = bin(boundary2d, RectangleBinner(500))


# function bin(boudary::Boundary{Dim, CRS, TP}, binner::AbstractBinner{Dim}) where {Dim, CRS, TP}
# hardclip: intersect(mesh, grid)