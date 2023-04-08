module UlamTypes

import MAT
import HDF5
import LinearAlgebra

export 
    UlamPolygon,
    UlamTrajectories,
    UlamDomain,
    UlamProblem,
    UlamInfo,
    UlamResult,
    PolyTable,
    InpolyResult,
    show

const global_poly_types::Vector{String} = ["sqr", "hex", "vor"]
const global_stoc_types::Vector{String} = ["data", "source"]
const global_traj_file_types::Vector{String} = ["mat", "h5"]
const global_poly_number_default::Int64 = 100

# sort edges if they aren't already

struct UlamPolygon 
    nodes::Matrix{Float64} 
    edges::Matrix{Int64} 
    center::Matrix{Float64} 
    poly_type::String

    function UlamPolygon(
        nodes::Matrix{<:Real};
        edges::Union{Matrix{<:Integer}, Nothing} = nothing, 
        center::Union{Matrix{<:Real}, Nothing} = nothing, 
        poly_type::String = "unk")

        @assert size(nodes, 1) > 0 
        @assert size(nodes, 2) == 2 
        @assert (poly_type in global_poly_types) || (poly_type == "unk")

        n_nodes = size(nodes, 1)

        if edges === nothing
            edges = [1:n_nodes;; [2:n_nodes; 1]] # assume provided nodes are sorted
        elseif edges != [1:n_nodes;; [2:n_nodes; 1]] # user provided unsorted nodes; sort them
            ss = sortslices(edges, dims = 1)
            order = [1]
            for i = 2:size(ss, 1)
                push!(order, ss[order[i - 1], 2])
            end

            nodes = nodes[order, :]
            edges = [1:n_nodes;; [2:n_nodes; 1]]
        end

        @assert size(edges) == size(nodes)

        if center === nothing
            center = [sum(nodes[:,1]) sum(nodes[:,2])]/n_nodes
        end

        @assert size(center) == (1, 2)

        new(nodes, edges, center, poly_type)
    end
end

function Base.show(io::IO, x::UlamPolygon)
    print(io, "UlamPolygon[")
    show(io, size(x.edges, 1))
    print("]")
end

function Base.show(io::IO, ::MIME"text/plain", x::Vector{UlamPolygon})
    print(io, "Vector{UlamPolygon}[")
    show(io, length(x))
    print("]")
end



struct UlamTrajectories 
    x0::Vector{Float64}
    y0::Vector{Float64}
    xT::Vector{Float64}
    yT::Vector{Float64}

    function UlamTrajectories(
        ;
        x0::Vector{<:Real},
        y0::Vector{<:Real}, 
        xT::Vector{<:Real}, 
        yT::Vector{<:Real})

        @assert length(x0) == length(y0) == length(xT) == length(yT) > 0
        new(x0, y0, xT, yT)
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

        @assert length(x0) == length(y0) == length(xT) == length(yT) > 0

        new(x0, y0, xT, yT)
    end
end  


function Base.show(io::IO, x::UlamTrajectories)
    print(io, "UlamTrajectories[")
    show(io, length(x.x0))
    print("]")
end


struct UlamDomain
    domain::Union{UlamPolygon, Nothing}
    corners::Vector{Float64}
    poly_type::String
    poly_number::Int64
    stoc_type::String
    stoc_polygon::Union{UlamPolygon,Nothing}
    rseed::Int64

    function UlamDomain(
        xmin::Real,
        xmax::Real,
        ymin::Real,
        ymax::Real;
        domain::Union{UlamPolygon, Nothing} = nothing,
        poly_type::String = global_poly_types[1],
        poly_number::Int64 = global_poly_number_default,
        stoc_type::String = global_stoc_types[1],
        stoc_polygon::Union{UlamPolygon, Nothing} = nothing,
        rseed::Int64 = 123)

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

        new(domain, corners, poly_type, poly_number, stoc_type, stoc_polygon), rseed
    end
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


struct UlamInfo
    n_polys_requested::Int64
    n_polys_no_data::Int64
    n_polys_dis::Int64
    n_counts_total::Int64
    n_counts_removed_scc::Int64
    poly_type::String
end


struct UlamResult
    P_closed::Matrix{Float64}
    P_open::SubArray{Float64, 2, Matrix{Float64}, Tuple{UnitRange{Int64}, UnitRange{Int64}}, false}
    pi_closed::Vector{Float64}
    pi_open::SubArray{Float64, 1, Vector{Float64}, Tuple{UnitRange{Int64}}, true}
    polys::Vector{UlamPolygon}
    polys_dis::Vector{UlamPolygon}
    counts::Vector{Int64}
    info::UlamInfo

    function UlamResult(
        P_closed::Matrix{Float64},
        polys::Vector{UlamPolygon},
        polys_dis::Vector{UlamPolygon},
        counts::Vector{Int64},
        info::UlamInfo)

        @assert size(P_closed, 1) == size(P_closed, 2)
        @assert size(P_closed, 1) > 0

        @assert length(polys) > 0
        
        pi_closed = abs.(LinearAlgebra.normalize(LinearAlgebra.eigvecs(transpose(P_closed))[:,size(P_closed)[1]], 1))

        P_open = view(P_closed, 1:size(P_closed, 1) - 1, 1:size(P_closed, 1) - 1)
        pi_open = view(pi_closed, 1:size(pi_closed, 1) - 1)

        new(P_closed, P_open, pi_closed, pi_open, polys, polys_dis, counts, info)
    end
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


struct PolyTable
    nodes::Matrix{Float64}
    edges::Matrix{Int64}
    n_polys::Int64

    function PolyTable(UPpolys::Vector{UlamPolygon})
        n_polys = length(UPpolys)
        @assert n_polys > 0

        n_nodes = sum(size(UPpolys[i].nodes, 1) for i = 1:n_polys)

        nodes = Matrix{Float64}(undef, n_nodes, 2)
        edges = Matrix{Int64}(undef, n_nodes, 3)

        counter_edge = 1
        counter_bot = 1

        for poly in UPpolys
            counter_top = counter_bot + size(poly.nodes, 1) - 1

            nodes[counter_bot:counter_top, :] = poly.nodes
            edges[counter_bot:counter_top, 1:2] = poly.edges .+ counter_bot .- 1
            edges[counter_bot:counter_top, 3] .= counter_edge

            counter_bot = counter_top + 1
            counter_edge = counter_edge + 1
        end

        new(nodes, edges, n_polys)
    end
end

function Base.show(io::IO, x::PolyTable)
    print(io, "PolyTable[")
    show(io, x.n_polys)
    print("]")
end

struct InpolyResult
    inds::Vector{Int64}
    contains::Dict{Int64, Int64}
end

end # module

