using LinearAlgebra

function SEBA(V, Rinit=nothing)

# Inputs: 
# V is pxr matrix (r vectors of length p as columns)
# Rinit is an (optional) initial rotation matrix.

# Outputs:
# S is pxr matrix with columns approximately spanning the column space of V
# R is the optimal rotation that acts on V, which followed by thresholding, produces S

    maxiter = 5000   #maximum number of iterations allowed
    F = qr(V) # Enforce orthonormality
    V = Matrix(F.Q)
    p, r = size(V)
    μ = 0.99 / sqrt(p)

    S = zeros(size(V))
    # Perturb near-constant vectors
    for j = 1:r
        if maximum(V[:,j]) - minimum(V[:,j]) < 1e-14
            V[:,j] = V[:,j] .+ (rand(p, 1) .- 1 / 2) * 1e-12
        end
    end

    # Initialise rotation
    if Rinit ≡ nothing
        Rnew = Matrix(I, r, r)
    else
        # Ensure orthonormality of Rinit
        F = svd(Rinit)
        Rnew = F.U * F.Vt
    end

    R = zeros(r, r)
    iter = 0
    while norm(Rnew - R) > 1e-14 && iter < maxiter
        iter = iter + 1
        R = Rnew
        Z = V * R'
        # Threshold to solve sparse approximation problem
        for i = 1:r
            Si = sign.(Z[:,i]) .* max.(abs.(Z[:,i]) .- μ, zeros(p))
            S[:,i] = Si / norm(Si)
        end
        # Polar decomposition to solve Procrustes problem
        F = svd(S' * V, full=false)
        Rnew = F.U * F.Vt
    end

    # Choose correct parity of vectors and scale so largest value is 1
    for i = 1:r
        S[:,i] = S[:,i] * sign(sum(S[:,i]))
        S[:,i] = S[:,i] / maximum(S[:,i])
    end

    # # Sort so that most reliable vectors appear first
    # ind = sortperm(vec(minimum(S, dims=1)), rev=true)
    # S = S[:, ind]

    return S, R

end

function attractors(P_closed, n_evec)
    evecs = eigen(transpose(P_closed)).vectors[:,end-n_evec+1:end] 
    # evecs = eigen(transpose(P_closed)).vectors[:,end-1-n_evec+1:end-1] # n_evec subdominant eigenfuctions

    if maximum(imag(evecs)) >= 10.0e-12
        # display(maximum(imag(evecs)))
        display("Warning: eigenvectors from SEBA are not real. Taking real part.")
    end

    S, R = SEBA(real(evecs))
    return S
end

function basins(P_closed, n_evec)
    evecs = eigen(P_closed).vectors[:,end-n_evec+1:end] 
    # evecs = eigen(P_closed).vectors[:,end-1-n_evec+1:end-1] # n_evec subdominant eigenfuctions

    if maximum(imag(evecs)) >= 10.0e-12
        # display(maximum(imag(evecs)))
        display("Warning: eigenvectors from SEBA are not real. Taking real part.")
    end

    S, R = SEBA(real(evecs))
    return S
end