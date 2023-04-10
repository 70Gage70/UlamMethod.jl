struct PolyTable{T<:Real, U<:Integer}
    nodes::Matrix{T}
    edges::Matrix{U}
    n_polys::U
end

function PolyTable(UPpolys::Vector{UlamPolygon{T}}) where {T<:Real}
    n_polys = length(UPpolys)
    if n_polys == 0 
        return PolyTable{T, Int64}(Matrix{T}(undef, 0, 2), Matrix{Int64}(undef, 0, 3), n_polys) 
    end

    n_nodes = sum(size(UPpolys[i].nodes, 1) for i = 1:n_polys)

    nodes = Matrix{T}(undef, n_nodes, 2)
    edges = Matrix{Int64}(undef, n_nodes, 3)

    counter_edge = 1
    counter_bot = 1

    for poly in UPpolys
        counter_top = counter_bot + size(poly.nodes, 1) - 1

        # The edges of an UlamPolygon are always sorted, so we know that the table of edges
        # is a matrix with size(poly.nodes, 1) rows and 2 columns. 
        poly_edges = [1:size(poly.nodes, 1);; [2:size(poly.nodes, 1); 1]]

        nodes[counter_bot:counter_top, :] = poly.nodes
        edges[counter_bot:counter_top, 1:2] = poly_edges .+ counter_bot .- 1
        edges[counter_bot:counter_top, 3] .= counter_edge

        counter_bot = counter_top + 1
        counter_edge = counter_edge + 1
    end

    return PolyTable{T, Int64}(nodes, edges, n_polys)
end

function Base.show(io::IO, x::PolyTable)
    print(io, "PolyTable[")
    show(io, x.n_polys)
    print("]")
end