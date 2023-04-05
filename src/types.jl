const global_bin_types = ["reg", "hex", "vor"]
const global_poly_types = ["reg", "hex", "vor", "unk"]

abstract type AbstractInPolygonCompatible end

struct UlamPolygon <: AbstractInPolygonCompatible
    nodes::Matrix{Float64} 
    edges::Matrix{Int64} 
    center::Matrix{Float64} 
    polys_type::String

    function UlamPolygon(
        nodes::Matrix{Float64};
        edges::Union{Matrix{Int64}, Nothing} = nothing, 
        center::Union{Matrix{Float64}, Nothing} = nothing, 
        polys_type::String = "unk")

        @assert size(nodes, 1) > 0 
        @assert size(nodes, 2) == 2 
        @assert polys_type in global_poly_types

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


struct PolyTable <: AbstractInPolygonCompatible
    nodes::Matrix{Float64}
    edges::Matrix{Int64}

    function PolyTable(UPpolys::Vector{UlamPolygon})
        @assert length(UPpolys) > 0

        n_nodes = sum(size(UPpolys[i].nodes, 1) for i = 1:length(UPpolys))

        nodes = Matrix{Float64}(undef, n_nodes, 2)
        edges = Matrix{Int64}(undef, n_nodes, 3)

        counter_edge = 1
        counter_bot = 1
        counter_top = counter_bot + size(poly.nodes, 1) - 1
        offset = size(UPpolys[1].nodes, 1)
        counter_top = counter_bot + offset - 1

        for i = 1:length(UPpolys)
            nodes[counter_bot:counter_top, :] = poly.nodes
            edges[counter_bot:counter_top, 1:2] = poly.edges .+ counter_top .- offset
            edges[counter_bot:counter_top, 3] .= counter_edge

            counter_bot = counter_top + 1
            counter_edge = counter_edge + 1
        end

        new(nodes, edges)
    end
end


struct UlamCovering
    polys::Vector{UlamPolygon}
    contains_data::Vector{Bool}
    contains_scc::Vector{Bool}
    contains_A::Vector{Bool}
    contains_B::Vector{Bool}


    function UlamCovering(
        polys::Vector{UlamPolygon};
        contains_data::Union{Vector{Bool}, Nothing} = nothing,
        contains_scc::Union{Vector{Bool}, Nothing} = nothing,
        contains_A::Union{Vector{Bool}, Nothing} = nothing,
        contains_B::Union{Vector{Bool}, Nothing}= nothing)

        @assert length(polys) > 0 

        if contains_data === nothing
            contains_data = fill(true, length(polys))
        end

        if contains_scc === nothing
            contains_scc = fill(true, length(polys))
        end

        if contains_A === nothing
            contains_A = fill(false, length(polys))
        end

        if contains_B === nothing
            contains_B = fill(false, length(polys))
        end

        new(polys, contains_data, contains_scc, contains_A, contains_B)
    end

end

include("binner-square.jl")
n_polys_test = 100
polys_test, polys_centers_test = square_binner(n_polys_test, [0, 1, 0, 1])
include("ulam-nirvana.jl")
iptest = inpoly_preprocess(polys_test)
uptest = [UlamPolygon(iptest["nodes"][i:i+3,:]) for i=1:4:4*n_polys_test]
uctest = UlamCovering(uptest)
pttest = PolyTable(uctest.polys)


##############################################################################################################################

# To define an Ulam method, need corners, x0y0xTyT data, binning type and number, stochasticization type 

struct UlamBoundary 
    xmin::Float64
    xmax::Float64
    ymin::Float64
    ymax::Float64

    function UlamBoundary(;xmin::R1, xmax::R2, ymin::R3, ymax::R4) where {R1 <: Real, R2 <: Real, R3 <: Real, R4 <: Real}
        @assert xmin < xmax
        @assert ymin < ymax
        new(xmin, xmax, ymin, ymax)
    end
end  

struct UlamTrajectories 
    x0::Vector{Float64}
    xT::Vector{Float64}
    y0::Vector{Float64}
    yT::Vector{Float64}

    function UlamTrajectories(;x0::Vector{R1}, y0::Vector{R2}, xT::Vector{R3}, yT::Vector{R4}) where {R1 <: Real, R2 <: Real, R3 <: Real, R4 <: Real}
        @assert length(x0) == length(y0) == length(xT) == length(yT) > 0
        new(x0, xT, y0, yT)
    end
end  

struct UlamBins 
    bin_type::String
    bin_number::Integer

    function UlamBins(;bin_type::String, bin_number::Integer)
        @assert bin_type in global_bin_types
        @assert bin_number > 0
        new(bin_type, bin_number)
    end
end  

struct UlamStocData
    stoc_type::String
    bin_number::Integer

    function UlamBins(;bin_type::String, bin_number::Integer)
        @assert bin_type in global_bin_types
        @assert bin_number > 0
        new(bin_type, bin_number)
    end
end  



struct UlamProblem
    corners::Vector{T} where {T<:Real}

    function UlamProblem(corners::Vector{T}) where {T<:Real}
        if length(corners) != 4
            throw(ArgumentError("Vector must have length 4"))
        else
            new(corners)
        end
    end
end