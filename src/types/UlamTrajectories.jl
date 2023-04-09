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