module DiscreteVectorBundles

# Write your package code here.

export LocalTriv
export CechCocycle

export bestorthtrans
export approxcocycledeath
export approxcechcocycle
export sw1
export eu
export orientbundle

export cechsimplex
export transitionfunction

include("ApproxCocycles.jl")


export chaintovector
export one∂
export zeroδ
export linearsolvemodp
export evaluate
include("FpLinearAlgebra.jl")

export geometricorientation
export μ
include("EigenBundle.jl")

export orientbundle!
export fliporientation!
function orientbundle!(localtriv::LocalTriv,threshold::Number)
    #take a local trivialisation of a discrete vector bundle and 
    #orient it at the time (threshold) in the filtration

    cocycle = approxcechcocycle(localtriv,threshold)
    w1 = chaintovector(localtriv.complex,sw1(cocycle),2)
    cobound = zeroδ(localtriv.complex,2,threshold)
    Θ = linearsolvemodp(cobound,w1,2)
    
    for i in eachindex(Θ)
        if Θ[i] == 1
            ind = index(localtriv.complex[0][i])
            b = deepcopy(localtriv.bases[ind])
            localtriv.bases[ind][:,2] = b[:,1]
            localtriv.bases[ind][:,1] = b[:,2]
        end
    end
end

function fliporientation!(localtriv::LocalTriv)
    #take a local trivialisation of a discrete vector bundle and invert its orientation
    for i in eachindex(localtriv.bases)
        ind = index(localtriv.complex[0][i])
        b = deepcopy(localtriv.bases[ind])
        localtriv.bases[ind][:,2] = b[:,1]
        localtriv.bases[ind][:,1] = b[:,2]
    end
end

#Placeholder for the full cohomology solver function
function solveincohomologybasis()
    return 0
end


end


