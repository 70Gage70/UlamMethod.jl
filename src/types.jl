module UlamTypes

import MAT
import HDF5
import LinearAlgebra

export 
    UlamPolygon,
    UlamTrajectories,
    UlamDomain,
    UlamInfo,
    UlamResult,
    PolyTable,
    InpolyResult,
    show

const global_poly_types::Vector{String} = ["sqr", "hex", "vor"]
const global_stoc_types::Vector{String} = ["data", "source"]
const global_traj_file_types::Vector{String} = ["mat", "h5"]
const global_poly_number_default::Dict = Dict("sqr" => 500, "hex" => 500, "vor" => 100)



struct UlamPolygon{T<:Real} 
    nodes::Matrix{T} 
    center::Matrix{T} 
end

function UlamPolygon(
    nodes::Matrix{<:Real};
    center::Union{Matrix{<:Real}, Nothing} = nothing, 
    edges::Union{Matrix{<:Integer}, Nothing} = nothing)

    @assert size(nodes, 1) > 0 
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


struct UlamTrajectories{T<:Real}
    x0::Vector{T}
    y0::Vector{T}
    xT::Vector{T}
    yT::Vector{T}
end  


function UlamTrajectories(
    ;
    x0::Vector{<:Real},
    y0::Vector{<:Real}, 
    xT::Vector{<:Real}, 
    yT::Vector{<:Real})

    @assert length(x0) == length(y0) == length(xT) == length(yT) > 0

    convert(Vector{Float64}, x0)
    promote(x0, y0, xT, yT)

    return UlamTrajectories{Float64}(x0, y0, xT, yT)
end


function UlamTrajectories(infile::String)
    extension = infile[findlast(==('.'), infile)+1:end]

    @assert extension in global_traj_file_types

    if extension == "mat"
        data_T = MAT.matopen(infile)
    elseif extension == "h5"
        data_T = HDF5.h5open(infile)
    end

    x0, y0, xT, yT = map(vec, read(data_T, "x0", "y0", "xT", "yT"))
    close(data_T)

    @assert eltype(x0) <: Real 
    @assert eltype(y0) <: Real 
    @assert eltype(xT) <: Real 
    @assert eltype(yT) <: Real 

    return UlamTrajectories(x0 = x0, y0 = y0, xT = xT, yT = yT)
end


function Base.show(io::IO, x::UlamTrajectories{T}) where {T<:Real}
    print(io, "UlamTrajectories[")
    show(io, length(x.x0))
    print("]")
end


struct UlamDomain{S<:AbstractString, T<:Real, U<:Integer}
    corners::Vector{T}
    domain::Union{UlamPolygon{T}, Nothing}
    poly_type::S
    poly_number::U
    stoc_type::S
    stoc_polygon::Union{UlamPolygon{T},Nothing}
    rseed::U
end

function UlamDomain(
    xmin::T,
    xmax::T,
    ymin::T,
    ymax::T;
    domain::Union{UlamPolygon{T}, Nothing} = nothing,
    poly_type::S = global_poly_types[1],
    poly_number::U = global_poly_number_default[poly_type],
    stoc_type::S = global_stoc_types[1],
    stoc_polygon::Union{UlamPolygon{T}, Nothing} = nothing,
    rseed::U = 123) where {S<:AbstractString, T<:Real, U<:Integer}

    @assert xmax > xmin
    @assert ymax > ymin

    corners::Vector{Float64} = [xmin, xmax, ymin, ymax]

    # if corners is ever replaced by a generic domain
    # if domain === nothing
    #     domain = UlamPolygon([xmin ymin; xmin ymax; xmax ymax; xmin ymax])
    # end

    @assert poly_type in global_poly_types
    @assert poly_number > 1
    @assert stoc_type in global_stoc_types

    if stoc_polygon != nothing && stoc_type != "source"
        # assume the user wants to use the source algorithm if they provide a stoc_polygon
        stoc_type = "source"
    end

    return UlamDomain{String, Float64, Int64}(
        corners,
        domain,
        poly_type,
        poly_number,
        stoc_type,
        stoc_polygon,
        rseed)
end

function Base.show(io::IO, x::UlamDomain)
    print(io, "UlamDomain[")
    show(io, x.corners)
    print(io, ", ")
    show(io, x.poly_number)
    print(io, "@")
    show(io, x.poly_type)
    if x.domain != nothing
        print(io, " w/ domain")
    end
    print("]")
end


struct UlamInfo{S<:AbstractString, U<:Integer}
    n_polys_requested::U
    n_polys_no_data::U
    n_polys_dis::U
    n_counts_total::U
    n_counts_removed_scc::U
    poly_type::S

    function UlamInfo(
        n_polys_requested::U,
        n_polys_no_data::U,
        n_polys_dis::U,
        n_counts_total::U,
        n_counts_removed_scc::U,
        poly_type::S) where {S<:AbstractString, U<:Integer}
    
        new{String, Int64}(
        n_polys_requested,
        n_polys_no_data,
        n_polys_dis,
        n_counts_total,
        n_counts_removed_scc,
        poly_type)
    end
end




struct UlamResult{S<:AbstractString, T<:Real, U<:Integer}
    P_closed::Matrix{T}
    P_open::SubArray{T, 2, Matrix{T}, Tuple{UnitRange{U}, UnitRange{U}}, false}
    pi_closed::Vector{T}
    pi_open::SubArray{T, 1, Vector{T}, Tuple{UnitRange{U}}, true}
    polys::Vector{UlamPolygon{T}}
    polys_dis::Vector{UlamPolygon{T}}
    counts::Vector{U}
    info::UlamInfo{S, U}
end


function UlamResult(
    P_closed::Matrix{T},
    polys::Vector{UlamPolygon{T}},
    polys_dis::Vector{UlamPolygon{T}},
    counts::Vector{U},
    info::UlamInfo{S, U}) where {S<:AbstractString, T<:Real, U<:Integer}

    @assert size(P_closed, 1) == size(P_closed, 2)
    @assert size(P_closed, 1) > 0

    @assert length(polys) > 0
    
    pi_closed = abs.(LinearAlgebra.normalize(LinearAlgebra.eigvecs(transpose(P_closed))[:,size(P_closed)[1]], 1))

    P_open = view(P_closed, 1:size(P_closed, 1) - 1, 1:size(P_closed, 1) - 1)
    pi_open = view(pi_closed, 1:size(pi_closed, 1) - 1)

    return UlamResult{S, T, U}(P_closed, P_open, pi_closed, pi_open, polys, polys_dis, counts, info)
end


function Base.show(io::IO, x::UlamResult)
    print(io, " UlamResult")
    println(io)
    print("  Polys requested: ")
    show(io, x.info.n_polys_requested)
    println(io)
    print("  Polys type: ")
    show(io, x.info.poly_type)
    println(io)
    print("  Polys obtained: ")
    show(io, length(x.polys))
    println(io)
    print("  Polys disconnected: ")
    show(io, length(x.polys_dis))
end


struct PolyTable{T<:Real, U<:Integer}
    nodes::Matrix{T}
    edges::Matrix{U}
    n_polys::U
end

function PolyTable(UPpolys::Vector{UlamPolygon{T}}) where {T<:Real}
    n_polys = length(UPpolys)
    @assert n_polys > 0

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

struct InpolyResult{U<:Integer}
    inds::Vector{U}
    contains::Dict{U, U}
end

end # module

