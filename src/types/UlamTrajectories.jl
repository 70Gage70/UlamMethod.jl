"""
    UlamTrajectories{T}

A container for trajectory data.
"""
struct UlamTrajectories{T<:Real}
    x0::Vector{T}
    y0::Vector{T}
    xT::Vector{T}
    yT::Vector{T}
end  

"""
    UlamTrajectories(;x0, y0, xT, yT)

Construct a container for trajectory data.
"""
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

"""
    UlamTrajectories(infile, [...])

Construct a container for trajectory data, loading it from `infile`, which should be a `.h5` or `.mat` file.

    ### Optional Arguments
    These are used if trajectory data in `infile` are named something other than `"x0"`, `"y0"`, `"xT"`, `"yT"`.
    - `x0_alias`
    - `y0_alias`
    - `xT_alias`
    - `yT_alias`
"""
function UlamTrajectories(
    infile::String; 
    x0_alias::String = "x0",
    y0_alias::String = "y0",
    xT_alias::String = "xT",
    yT_alias::String = "yT")

    extension = infile[findlast(==('.'), infile)+1:end]

    @assert extension in global_traj_file_types "Require a .mat or .h5 file."

    if extension == "mat"
        data_T = MAT.matopen(infile)
    elseif extension == "h5"
        data_T = HDF5.h5open(infile)
    end

    @assert x0_alias in collect(keys(data_T)) "$(x0_alias) not in $(infile)"
    @assert y0_alias in collect(keys(data_T)) "$(y0_alias) not in $(infile)"
    @assert xT_alias in collect(keys(data_T)) "$(xT_alias) not in $(infile)"
    @assert yT_alias in collect(keys(data_T)) "$(yT_alias) not in $(infile)"

    x0, y0, xT, yT = map(vec, read(data_T, x0_alias, y0_alias, xT_alias, yT_alias))
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