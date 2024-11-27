"""
    struct Boundary{Dim, M, CRS}

A container type for the computational domain within which Ulam's method is applied.

The boundary has dimension `Dim` embedded in manifold `M` with coordinate reference system `CRS`.

The boundary is partitioned according to an [`BinningAlgorithm`](@ref).

### Fields

- `boundary`: A `Polytope` defining the boundary.

### Methods

    points(boundary)

Calculate a `Dim x N` matrix of boundary vertices.

### 1D Constructors

    Boundary(x_start, x_end)

The boundary is a continuous line segment `Segment` between `x_start` and `x_end`.

### 2D Constructors

    Boundary(verts)

The boundary is a closed polygon `Ngon` with vertices `verts`, which should be be a
vector of `(x, y)` coordinates or a `2 x N` matrix.

    Boundary(xmin, xmax, ymin, ymax)

A convenience constructor is also provided for the case of a rectangular boundary.

### ≥3D Constructors

    Boundary(corner_min, corner_max)

The boundary is a `HyperRectangle` with minimum and maximum vertices `corner_min`, `corner_max`
where the corners are `NTuple`s, `(x_min, y_min, z_min, ...)`, `(x_max, y_max, z_max, ...)`.

### Automatric Boundary Construction

    AutoBoundary(traj; nirvana = 0.10)

Construct a rectangular boundary in arbitrary dimensions based on the [`Trajectories`](@ref) in `traj`.

The boundary is placed such that a fraction `nirvana` of the data is in nirvana. For example, with the default 
value of `0.10`, roughly `10%` of all of the datapoints (split between `x0` and `xT`) will be outside the boundary with 
a roughly equal amount on each side.

The shape of the boundary can be further controlled by providing `nirvana` as a vector of length `Dim` of tuples of the form `(min, max)` such 
that a fraction `min` (respectively `max`) will be "below" (respectively "above") the boundary along each dimension. 

    AutoBoundary2D(traj; nirvana = 0.10, tol = 0.001)

Construct a tight boundary in two dimensions based on the [`Trajectories`](@ref) in `traj`.

The boundary in question is formed by scaling the convex hull of the trajectory data until a fraction `nirvana` of 
the data are bouside the boundary to within `tol` percent. 
"""
struct Boundary{Dim, M, CRS}
    boundary::Polytope{Dim, M, CRS}

    function Boundary(; boundary::Polytope{Dim, M, CRS}) where {Dim, M, CRS}
        return new{Dim, M, CRS}(boundary)
    end
end

### 1D
function Boundary(x_start::Real, x_end::Real)
    @argcheck x_start < x_end

    return Boundary(boundary = Segment(Point(float(x_start)), Point(float(x_end))))
end

function points(boundary::Boundary{1, M, CRS}) where {M, CRS}
    verts = boundary.boundary.vertices
    
    return stack([[coords(v).x.val] for v in verts])
end

### 2D
function Boundary(verts::Vector{<:Tuple{Real, Real}})
    @argcheck all(length.(verts) .== 2)

    return Boundary(boundary = Ngon(verts...))
end

function Boundary(verts::Matrix{<:Real})
    @argcheck size(verts, 1) == 2

    return Boundary([(verts[1, i], verts[2, i]) for i = 1:size(verts, 2)])
end

function Boundary(xmin::Real, xmax::Real, ymin::Real, ymax::Real)
    @argcheck xmin < xmax
    @argcheck ymin < ymax

    return Boundary([(xmin, ymin), (xmax, ymin), (xmax, ymax), (xmin, ymax)])
end

function points(boundary::Boundary{2, M, CRS}) where {M, CRS}
    verts = boundary.boundary.vertices
    
    return stack([[coords(v).x.val, coords(v).y.val] for v in verts])
end

### ≥3D
function Boundary(corner_min::NTuple{N, Real}, corner_max::NTuple{N, Real}) where {N}
    @argcheck length(corner_min) == length(corner_max)
    @argcheck N ≥ 3 "Prefer to not use a HyperRectangle in 1 or 2 dimensions"
    return Boundary(boundary = HyperRectangle(corner_min, corner_max))
end

### AutoBoundary
function AutoBoundary(
    traj::Trajectories{Dim}; 
    nirvana::Union{Real, Vector{<:Tuple{Real, Real}}} = 0.10) where {Dim}

    if nirvana isa Real
        @argcheck 0 <= nirvana < 1
        nirvana = [(nirvana/(2*Dim), nirvana/(2*Dim)) for _ = 1:Dim]
    else
        @argcheck length(nirvana) == Dim
        facs = nirvana |> Iterators.flatten |> collect
        @argcheck all(0 .<= facs .<= 1)
        @argcheck sum(facs) <= 1
    end

    data = [traj.x0 ;; traj.xT]
    n_points = [size(data, 2) .* fac .|> x -> max(1, ceil(Int64, x)) for fac in nirvana]
    
    sort_perms = [sortperm(data[d,:]) for d = 1:Dim]
    bounds = [(
        data[d, sort_perms[d][n_points[d][1]]], 
        data[d, sort_perms[d][end-n_points[d][2]]]) for d = 1:Dim]

    if Dim == 1
        return Boundary(bounds[1][1], bounds[1][2])
    elseif Dim == 2
        return Boundary(bounds[1][1], bounds[1][2], bounds[2][1], bounds[2][2])
    else
        corner_min = [bounds[i][1] for i = 1:Dim] |> Tuple
        corner_max = [bounds[i][2] for i = 1:Dim] |> Tuple
        return Boundary(corner_min, corner_max)
    end
end

function AutoBoundary2D(
    traj::Trajectories{Dim}; 
    nirvana::Real = 0.10,
    tol::Real = 0.001) where {Dim}

    @argcheck Dim == 2
    @argcheck 0 <= nirvana <= 1
    @argcheck tol > 0

    # compute convex hull of data
    data = [traj.x0 ;; traj.xT]
    ps = [Point(data[:,i]...) for i = 1:size(data, 2)] |> PointSet
    hull = convexhull(ps).rings[1].vertices |> x -> Ngon(x...)
    hull = Scale(1.01)(hull) # scale up slightly to avoid floating point errors at the boundary

    # compute the fraction of data inside the hull when it has been scaled by `scale`
    _interior(scale) = _membership2d(data, Bins([Scale(scale)(hull)])) |> x -> length(findall(!isnothing, x))/size(data, 2)
    
    interior = 1 - nirvana
    a, b = 0.0, 1.0
    c = (a + b)/2
    n_iter = 0
    
    # use bisection to find the correct scale
    while abs(_interior(c) - interior)/interior > tol
        if _interior(c) - interior > 0
            b = c
        else
            a = c
        end
        
        c = (a + b)/2
        n_iter += 1
    
        if n_iter >= 20
            break
        end
    end

    return Boundary(boundary = Scale(c)(hull))
end