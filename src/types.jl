module UlamTypes

import MAT
import HDF5
import LinearAlgebra

export 
    UlamPolygon,
    PolyTable,
    # UlamCovering,
    UlamTrajectories,
    UlamDomain,
    UlamProblem

const global_poly_types::Vector{String} = ["reg", "hex", "vor"]
const global_stoc_types::Vector{String} = ["data", "source"]
const global_traj_file_types::Vector{String} = ["mat", "h5"]
const global_poly_number_default::Int64 = 100


struct UlamPolygon 
    nodes::Matrix{Float64} 
    edges::Matrix{Int64} 
    center::Matrix{Float64} 
    polys_type::String

    function UlamPolygon(
        nodes::Matrix{<:Real};
        edges::Union{Matrix{<:Integer}, Nothing} = nothing, 
        center::Union{Matrix{<:Real}, Nothing} = nothing, 
        polys_type::String = "unk")

        @assert size(nodes, 1) > 0 
        @assert size(nodes, 2) == 2 
        @assert (polys_type in global_poly_types) || (polys_types == "unk")

        n_nodes = size(nodes, 1)

        if edges === nothing
            edges = [1:n_nodes;; [2:n_nodes; 1]]
        end

        @assert size(edges) == size(nodes)

        if center === nothing
            center = [sum(nodes[:,1]) sum(nodes[:,2])]/n_nodes
        end

        @assert size(center) == (1, 2)

        new(nodes, edges, center, polys_type)
    end
end


struct PolyTable
    nodes::Matrix{Float64}
    edges::Matrix{Int64}

    function PolyTable(UPpolys::Vector{UlamPolygon})
        @assert length(UPpolys) > 0

        n_nodes = sum(size(UPpolys[i].nodes, 1) for i = 1:length(UPpolys))

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

        new(nodes, edges)
    end
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


struct UlamDomain
    domain::UlamPolygon
    corners::Vector{Float64}
    poly_type::String
    poly_number::Int64
    stoc_type::String
    stoc_polygon::Union{UlamPolygon,Nothing}

    function UlamDomain(
        xmin::Real,
        xmax::Real,
        ymin::Real,
        ymax::Real;
        domain::Union{UlamPolygon, Nothing} = nothing,
        poly_type::String = global_poly_types[1],
        poly_number::Int64 = global_poly_number_default,
        stoc_type::String = global_stoc_types[1],
        stoc_polygon::Union{UlamPolygon, Nothing} = nothing)

        @assert xmax > xmin
        @assert ymax > ymin

        corners::Vector{Float64} = [xmin, xmax, ymin, ymax]

        if domain === nothing
            domain = UlamPolygon([xmin ymin; xmin ymax; xmax ymax; xmin ymax])
        end

        @assert poly_type in global_poly_types
        @assert poly_number > 1
        @assert stoc_type in global_stoc_types

        new(domain, corners, poly_type, poly_number, stoc_type, stoc_polygon)
    end
end

# returndict = begin Dict(
#     "P_closed" => P_closed, 
#     "P_open" => P_open,
#     "pi_closed" => pi_closed,
#     "pi_open" => pi_open,
#     "counts" => final_counts,
#     "leaves" => leaves,
#     "polys" => vcells,
#     "polys_centers" => vcenters,
#     "polys_dis" => vcells_dis,
#     "info" => info)
# end

struct UlamCovering
    polys::Vector{UlamPolygon}
    counts::Vector{Int64}
    contains_scc::Vector{Bool}


    function UlamCovering(
        polys::Vector{UlamPolygon};
        counts::Vector{Int64} = fill(0, length(polys)),
        contains_scc::Vector{Bool} = fill(true, length(polys)))

        @assert length(polys) > 0 

        new(polys, counts, contains_scc)
    end

end

struct UlamResult
    covering::UlamCovering
    P_closed::Matrix{Float64}
    P_open::SubArray{Float64, 2, Matrix{Float64}, Tuple{UnitRange{Int64}, UnitRange{Int64}}, false}
    pi_closed::Vector{Float64}
    pi_open::SubArray{Float64, 1, Vector{Float64}, Tuple{UnitRange{Int64}}, true}
    info::String

    function UlamResult(
        covering::UlamCovering,
        P_closed::Matrix{Float64},
        info::String)

        @assert size(P_closed, 1) == size(P_closed, 2)
        @assert size(P_closed, 1) > 0
        
        pi_closed = abs.(normalize(LinearAlgebra.eigvecs(transpose(P_closed))[:,size(P_closed)[1]], 1))

        P_open = view(P_closed, 1:size(P_closed, 1) - 1, 1:size(P_closed, 1) - 1)
        pi_open = view(pi_open, 1:size(pi_open) - 1, 1)

        new(covering, P_closed, P_open, pi_closed, pi_open, info)
    end
end

end # module

# include("binner-square.jl")
# n_polys_test = 100
# polys_test, polys_centers_test = square_binner(n_polys_test, [0, 1, 0, 1])
# include("ulam-nirvana.jl")
# iptest = inpoly_preprocess(polys_test)
# uptest = [UlamPolygon(iptest["nodes"][i:i+3,:]) for i=1:4:4*n_polys_test]
# uctest = UlamCovering(uptest)
# pttest = PolyTable(uctest.polys)


##############################################################################################################################

# To define an Ulam method, need corners, x0y0xTyT data, binning type and number, stochasticization type 

# struct UlamBoundary 
#     xmin::Float64
#     xmax::Float64
#     ymin::Float64
#     ymax::Float64

#     function UlamBoundary(;xmin::R1, xmax::R2, ymin::R3, ymax::R4) where {R1 <: Real, R2 <: Real, R3 <: Real, R4 <: Real}
#         @assert xmin < xmax
#         @assert ymin < ymax
#         new(xmin, xmax, ymin, ymax)
#     end
# end  

# struct UlamTrajectories 
#     x0::Vector{Float64}
#     xT::Vector{Float64}
#     y0::Vector{Float64}
#     yT::Vector{Float64}

#     function UlamTrajectories(;x0::Vector{R1}, y0::Vector{R2}, xT::Vector{R3}, yT::Vector{R4}) where {R1 <: Real, R2 <: Real, R3 <: Real, R4 <: Real}
#         @assert length(x0) == length(y0) == length(xT) == length(yT) > 0
#         new(x0, xT, y0, yT)
#     end
# end  

# struct UlamBins 
#     bin_type::String
#     bin_number::Integer

#     function UlamBins(;bin_type::String, bin_number::Integer)
#         @assert bin_type in ["reg", "hex", "vor"]
#         @assert bin_number > 0
#         new(bin_type, bin_number)
#     end
# end  

# struct UlamStocData
#     stoc_type::String
#     bin_number::Integer

#     function UlamBins(;bin_type::String, bin_number::Integer)
#         @assert bin_type in ["reg", "hex", "vor"]
#         @assert bin_number > 0
#         new(bin_type, bin_number)
#     end
# end  



# struct UlamProblem
#     corners::Vector{T} where {T<:Real}

#     function UlamProblem(corners::Vector{T}) where {T<:Real}
#         if length(corners) != 4
#             throw(ArgumentError("Vector must have length 4"))
#         else
#             new(corners)
#         end
#     end
# end