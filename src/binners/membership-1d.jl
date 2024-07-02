function _membership1d(data::Matrix{<:Real}, bins::Bins{1, CRS}) where {CRS}
    @argcheck size(data, 1) == 1 
    bins = bins.bins
    memb = Union{Nothing, Int64}[]

    for x0 in data
        memb_i = nothing
        for i = 1:length(bins)
            bin_ = bins[i].vertices .|> x -> coords(x).x.val
            if bin_[1] < x0 < bin_[2]
                memb_i = i
                continue
            end
        end
        push!(memb, memb_i)
    end

    return memb
end

function _membership1d(traj::Trajectories{1}, bins::Bins{1, CRS}) where {CRS}  
    return (_membership1d(traj.x0, bins), _membership1d(traj.xT, bins))
end