struct Trajectories{Dim}
    x0::Matrix{Float64}
    xT::Matrix{Float64}

    function Trajectories(x0::Matrix{Float64}, xT::Matrix{Float64})
        @argcheck size(x0) == size(xT)
        @argcheck all(size(x0) .> 0)
        @argcheck size(x0, 1) in [1, 2] "Dimension of data must be 1 or 2."

        return new{size(x0, 1)}(x0, xT)
    end
end
