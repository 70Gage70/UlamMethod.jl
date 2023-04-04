# using Parameters

abstract type AbstractUlamCorners end
abstract type AbstractTrajectoryData end
abstract type AbstractBinData end

# To define an Ulam method, need corners, x0y0xTyT data, binning type and number, stochasticization type 

struct UlamCorners <: AbstractUlamCorners
    xmin::Float64
    xmax::Float64
    ymin::Float64
    ymax::Float64

    function UlamCorners(;xmin::Real, xmax::Real, ymin::Real, ymax::Real)
        @assert xmin < xmax
        @assert ymin < ymax
        new(xmin, xmax, ymin, ymax)
    end
end  


struct UlamTrajectories <: AbstractTrajectoryData
    x0::Vector{Float64}
    xT::Vector{Float64}
    y0::Vector{Float64}
    yT::Vector{Float64}

    function UlamTrajectories(;x0::Vector{Real}, y0::Vector{Real}, xT::Vector{Real}, yT::Vector{Real})
        @assert length(x0) == length(y0) == length(xT) == length(yT) > 0
        new(x0, xT, y0, yT)
    end
end  

struct UlamBins <: AbstractBinData
    bin_type::String
    bin_number::Integer

    function UlamBins(;bin_type::String, bin_number::Integer)
        @assert bin_type in ["reg", "hex", "vor"]
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