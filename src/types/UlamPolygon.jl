
struct UlamPolygon{T<:Real} 
    nodes::Matrix{T} 
    center::Matrix{T} 
end

function UlamPolygon(
    nodes::Matrix{<:Real};
    center::Union{Matrix{<:Real}, Nothing} = nothing, 
    edges::Union{Matrix{<:Integer}, Nothing} = nothing)

    @assert size(nodes, 1) > 2 "A polygon must have at least three nodes." 
    @assert size(nodes, 2) == 2 

    n_nodes = size(nodes, 1)

    if (edges != nothing) && (edges != [1:n_nodes;; [2:n_nodes; 1]]) # user provided unsorted nodes; sort them
        @assert size(edges) == size(nodes)

        ss = sortslices(edges, dims = 1)
        order = [1]
        for i = 2:size(ss, 1)
            push!(order, ss[order[i - 1], 2])
        end

        nodes = nodes[order, :]
    end

    if center === nothing
        center = [sum(nodes[:,1]) sum(nodes[:,2])]/n_nodes
    end

    @assert size(center) == (1, 2)

    nodes = convert(Matrix{Float64}, nodes)
    nodes, center = promote(nodes, center)

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