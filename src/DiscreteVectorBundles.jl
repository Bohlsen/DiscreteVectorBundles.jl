module DiscreteVectorBundles

using Nemo

#function exports
export LocalTriv
export filtration
export bases

export CechCocycle
export cechsimplex
export transitionfunction

export bestorthtrans
export approxcocycledeath
export approxcechcocycle
export sw1
export eu
export orientbundle
include("ApproxCocycles.jl")


export chaintovector
export one∂
export two∂
export zeroδ
export oneδ
export linearsolvemodp
export evaluate
include("FpLinearAlgebra.jl")

export geometricorientation
export μ
export constructrealeigenbundle
export constructcomplexeigenbundle
export canonicallyorientatedplane
include("EigenBundle.jl")

export orientbundle!
export fliporientation!
export solveincohomologybasis
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

function solveincohomologybasis(complex::Ripserer.AbstractFiltration,
                                cochain::Ripserer.Chain{<:Integer,Ripserer.Simplex{1,Float64,Int64}},
                                p::Number,
                                threshold::Number)
        #Take a 1 cochain and solve for its expression in the basis of cohomology produced by Ripserer's 
        #standard representative algorithm
        field, = residue_ring(ZZ,p)
        diagram = ripserer(complex,modulus = p,reps=true)
        H1basis = [c for c in reverse(diagram[2]) if death(c) > threshold]

        cocyclebasisvectors = matrix(field,zero_matrix(ZZ,length(complex[1]),0))
        for c in H1basis
            cocyclebasisvectors = matrix(field,[cocyclebasisvectors chaintovector(complex,representative(c),p,threshold)])
        end

        δ = zeroδ(complex,p,threshold)
        cocyclebasisvectors = matrix(field,[cocyclebasisvectors δ])
        b = chaintovector(complex,cochain,p,threshold)

        a = linearsolvemodp(cocyclebasisvectors,b,p)
        return a,H1basis
end

function solveincohomologybasis(complex::Ripserer.AbstractFiltration,
                                cochain::Ripserer.Chain{<:Integer,Ripserer.Simplex{2,Float64,Int64}},
                                p::Number,
                                threshold::Number)
        #Take a 2 cochain and solve for its expression in the basis of cohomology produced by Ripserer's 
        #standard representative algorithm

        field, = residue_ring(ZZ,p)
        diagram = ripserer(complex,dim_max=2,modulus = p,reps=true)
        H2basis = [c for c in reverse(diagram[3]) if death(c) > threshold]

        cocyclebasisvectors = matrix(field,zero_matrix(ZZ,length(complex[2]),0))
        for c in H2basis
            cocyclebasisvectors = matrix(field,[cocyclebasisvectors chaintovector(complex,representative(c),p,threshold)])
        end

        δ = oneδ(complex,p,threshold)
        cocyclebasisvectors = matrix(field,[cocyclebasisvectors δ])
        b = chaintovector(complex,cochain,p,threshold)

        a = linearsolvemodp(cocyclebasisvectors,b,p)
        return a,H2basis
end

end


