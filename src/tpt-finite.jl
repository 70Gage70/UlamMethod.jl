"""
Finite time TPT

Written by Gage Bonner September 2022
"""

###################################################################################################
###################################################################################################
###################################################################################################

using LinearAlgebra

"""
Subtracts x from y, setwise. 
e.g.  mysetdiff([1,2,3,4,5], [1,2,3]) = [4,5]
"""

function mysetdiff(y, x)
    res = Vector{eltype(y)}(undef, length(y) - length(x))
    i = 1
    @inbounds for el in y
        el ∈ x && continue
        res[i] = el
        i += 1
    end

    res
end

"""
Left eigenvector of P, normalized appropriately. The stationary distribution of P if it's row stochastic.
"""

function Ppi(P)
    return transpose(abs.(normalize(eigvecs(transpose(P))[:,size(P)[1]], 1)))
end

"""
In the case of the time homogenous-chain, this populates a 3d tensor 
where each "slice" is P. Note that [a ;;; b] constructs a 3d tensor
explicitly. The length of the time horizon is N.
"""

function constantP(N, P) 
    nStates = size(P)[1]
    q = zeros(nStates, nStates, N)
    for i = 1:N
        q[:, :, i] = P
    end
    return q
end

"""
The density of the chain at all times, obtained by repeated multiplication of the
transition matrix P by the initial density lambda0. Each row corresponds to a different
state in the chain and each column corresponds to increasing time.
"""

# rows are states, columns are n
function lambdan(lambda0, N, P)
    nStates = size(P)[1]
    L = zeros(nStates, N)
    L[:, 1] = lambda0

    for i=1:N - 1
        L[:, i + 1] = transpose(L[:, i])*P[:, :, i]
    end
    
    return L
end


"""
Forward committor. Solved iteratively with the boundary conditions that q+(B, N) = 1
and q+(A, N) = 0. The length of the time horizon is N. 

Note that lambda0 isn't actually used in qForward.
"""

function qForward(N, lambda0, ALLindices, Aindices, Bindices, P)
    Cindices = mysetdiff(ALLindices, [Aindices ; Bindices])
    nStates = size(P)[1]
    qForward = zeros(nStates,N)

    for i in Bindices # includes boundary condition
        for j in 1:N
            qForward[i, j] = 1 
        end
    end

    for n in N - 1:-1:1 # descending order since boundary condition is at n = N - 1
        for i in Cindices
            qForward[i, n] = sum([P[i, j, n]*qForward[j, n + 1] for j in ALLindices])
        end
    end
    
    return qForward
end

"""
Backward committor. Solved iteratively with the boundary conditions that q-(A, N) = 1
and q-(B, N) = 0. The length of the time horizon is N. This function technically
uses Pminus (the transition matrix of the backwards chain) but we don't need to 
compute it since each element can just be calculated when it's used in the iteration.

Note that lambda0 isn't actually used in qForward.
"""

function qBackward(N, lambda0, ALLindices, Aindices, Bindices, P)
    Cindices = mysetdiff(ALLindices, [Aindices ; Bindices])
    nStates = size(P)[1]
    qBackward = zeros(nStates,N)
    lambdanVALS = lambdan(lambda0, N, P)

    for i in Aindices # includes boundary condition
        for j in 1:N
            qBackward[i, j] = 1 
        end
    end

    for n in 2:1:N # for all times
        for i in Cindices # for all committors
            qBackward[i, n] = 0.0
            for j in ALLindices
                if lambdanVALS[i, n] > 0
                    qBackward[i, n] = qBackward[i, n] + (lambdanVALS[j, n - 1]/lambdanVALS[i, n])*P[j, i, n - 1]*qBackward[j, n - 1] 
                end
            end
        end
    end
    
    return qBackward
end

"""
All the statistics for finite time tpt, returned as a dict.
    "piStat" => stationary distribution [states]
    "q+" => forward committor [states, time]
    "q-" => backward committor [states, time] 
    "muAB" => un-normalized reactive density [states, time] 
    "muABnorm" => normalized reactive density [states, time] 
    "fij" => reactive current [state i, state j, time] 
    "f+" => effective reactive current [state i, state j, time]
    "kAB" => reactive rate [scalar] 
    "tAB" => reactive time [scalar]
"""

function tpt_finite_stats(N, lambda0, ALLindices, Aindices, Bindices, P)
    Cindices = mysetdiff(ALLindices, [Aindices ; Bindices])
    nStates = size(P)[1]
    lambdanVALS = lambdan(lambda0, N, P)
    qplus = qForward(N, lambda0, ALLindices, Aindices, Bindices, P)
    qminus = qBackward(N, lambda0, ALLindices, Aindices, Bindices, P)
    

    # reactive density
    muAB = zeros(nStates, N)
    
    for n = 2:N-1 # exclude places with mu = 0
        muAB[:, n] = qminus[:, n].*lambdanVALS[:, n].*qplus[:, n]
    end

    # normalization factor
    ZAB = [sum(muAB[:, n]) for n = 1:N]

    # normalized reactive density
    muABnorm = zeros(nStates, N)
    for n = 2:N-1 # exclude places with mu = 0
        muABnorm[:, n] = muAB[:, n]/ZAB[n]
    end    
    
    # reactive current
    fijAB = zeros(nStates, nStates, N)
    
    for n = 1:N-1 # exclude places with fij = 0
        for i in ALLindices
            for j in ALLindices
                fijAB[i, j, n] = qminus[i, n]*lambdanVALS[i, n]*P[i,j,n]*qplus[j, n+1]
            end
        end
    end
    
    # effective reactive current
    fijplus = zeros(nStates, nStates, N)
    
    for n = 1:N-1 # exclude places with fij = 0
        for i in ALLindices
            for j in ALLindices
                fijplus[i, j, n] = max(fijAB[i, j, n] - fijAB[j, i, n], 0)
            end
        end
    end
    
    # rate of leaving A
    kAout = zeros(N - 1)
    
    for n = 1:N - 1
        kAout[n] = sum(fijAB[i,j,n] for i in Aindices for j in ALLindices)
    end
    
    # overall rate
    kAB = (1/N)*sum(kAout)

    # overall time
    tAB = (1/N)*sum(ZAB)/kAB  

    tptDict = begin Dict(
        "piStat" => Ppi(P[:, :, 1]), 
        "q+" => qplus, 
        "q-" => qminus, 
        "muAB" => muAB, 
        "muABnorm" => muABnorm, 
        "fij" => fijAB, 
        "f+" => fijplus, 
        "kAB" => kAB, 
        "tAB" => tAB) 
    end

    return tptDict
           

end