"""
    UlamPolygon{T}

A polygon defined by a series of `nodes` and carrying a `center`.
"""
struct UlamPolygon{T<:Real} 
    nodes::Matrix{T} 
    center::Matrix{T} 
end

"""
    UlamPolygon(nodes; edges = nothing)

Construct an `UlamPolygon` based on a `nodes` matrix with `n` rows and `2` columns.

The `nodes` should not have a closing node, i.e. the last node should not be a repeat of the first.

The edges of the polygon defined by `nodes` are connected in order.

### Optional Arguments
- `edges`: an `n` by `2` matrix which specifies the edge connections betwen `nodes`. Used if nodes are not already sorted.
"""
function UlamPolygon(
    nodes::Matrix{<:Real};
    edges::Union{Matrix{<:Integer}, Nothing} = nothing)

    @assert size(nodes, 1) > 2 "A polygon must have at least three nodes." 
    @assert size(nodes, 2) == 2 

    n_nodes = size(nodes, 1)

    if (edges !== nothing) && (edges != [1:n_nodes;; [2:n_nodes; 1]]) # user provided unsorted nodes; sort them
        @assert size(edges) == size(nodes)

        ss = sortslices(edges, dims = 1)
        order = [1]
        for i = 2:size(ss, 1)
            push!(order, ss[order[i - 1], 2])
        end

        nodes = nodes[order, :]
    end

    nodes = convert(Matrix{Float64}, nodes)
    center = [sum(nodes[:,1]) sum(nodes[:,2])]/n_nodes

    return UlamPolygon{Float64}(nodes, center)
end

function Base.show(io::IO, x::UlamPolygon)
    print(io, "UlamPolygon[")
    show(io, size(x.nodes, 1))
    print("]")
end

function Base.show(io::IO, ::MIME"text/plain", x::Vector{UlamPolygon{T}}) where {T<:Real}
    print(io, "Vector{UlamPolygon}[")
    show(io, length(x))
    print("]")
end