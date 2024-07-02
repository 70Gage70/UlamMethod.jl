function Base.show(io::IO, x::Trajectories)
    dim = size(x.x0, 1)
    n_traj = size(x.x0, 2)
    print(io, "Trajectories{$(dim)}[$(n_traj)]")
end

function Base.show(io::IO, x::Boundary)
    dim = typeof(x).parameters[1]
    println(io, "Boundary{$(dim)}")
    println(io, "|")
    display(x.boundary)
end

function Base.show(io::IO, x::Bins)
    dim = typeof(x).parameters[1]
    n_bins = length(x.bins)
    print(io, "Bins{$(dim)}[$(n_bins)]")
end

function Base.show(io::IO, x::BinningAlgorithm)
    name_ba = string(typeof(x).name.name)
    name_bn = string(typeof(x.boundary.boundary).name.name)
    print(io, "$(name_ba)[$(name_bn) boundary, $(length(x.bins.bins)) bins]")
end

function Base.show(io::IO, x::ReinjectionAlgorithm)
    name_ba = string(typeof(x).name.name)
    print(io, "$(name_ba)")
end

function Base.show(io::IO, x::UlamResult)
    name_ba = string(typeof(x).name.name)
    dim = typeof(x).parameters[1]
    println(io, "$(name_ba){$(dim)D}")
    println(io, "├─$(x.binner)")
    print(io, "└─$(length(x.bins_dis.bins)) disconnected bins")
end
