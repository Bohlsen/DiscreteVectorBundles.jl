using Ripserer
using LinearAlgebra

struct LocalTriv
    #store both a full filtration over the point cloud and a local basis for the fiber at each vertex
    complex::Ripserer.AbstractFiltration{<:Integer,<:Number}
    bases::AbstractVector{<:AbstractVecOrMat}
end

struct CechCocycle
    simplex::Ripserer.Simplex
    Ω::Matrix{<:Number}
end
cechsimplex(elem::CechCocycle) = elem.simplex
transitionfunction(elem::CechCocycle) = elem.Ω


function bestorthtrans(s1::AbstractVecOrMat{<:Number},
                       s2::AbstractVecOrMat{<:Number})
        #Take two stiefel matrices (columns frame a vector space) and compute the 
        #smallest orthogonal matrix connecting them by solving Procustes problem

        u,s,v = svd(s1'*s2)
        return u*v'
end

function approxcocycledeath(
    localtriv::LocalTriv,
    tolerance::Number = 0.5,
    max_death::Number = Inf)
    #compute the maximum death time at which the local triviaisation generates a 
    #approximate cocycle with a given tolerance

    twoskeleton = sort(localtriv.complex[2]) #pull all 2 simplexes in the complex and sort by when they are added in the filtration

    for simplex in twoskeleton
        b = birth(simplex)
        if b > max_death
            return max_death
        else
            i,j,k = vertices(simplex)
            Ωij = bestorthtrans(localtriv.bases[i],localtriv.bases[j])
            Ωjk = bestorthtrans(localtriv.bases[j],localtriv.bases[k])
            Ωik = bestorthtrans(localtriv.bases[i],localtriv.bases[k])
            if norm(Ωik-Ωij*Ωjk) >= tolerance
                return b
            end
        end
    end
    return max_death
end


function approxcechcocyle(
    localtriv::LocalTriv,
    threshold::Number)
    #compute the approximate Cech cocycle from a local trivialisation
    #using simplices up to death=threshold

    cocycle = Vector{CechCocycle}(undef,0)
    oneskeleton = sort(localtriv.complex[1]) #pull all 1 simplexes in the complex and sort by when they are added in the filtration

    for edge in oneskeleton
        b = birth(edge)
        if b > threshold
            return cocycle
        end
        i,j = vertices(edge)
        Ωij=bestorthtrans(localtriv.bases[i],localtriv.bases[j])
        push!(cocycle,CechCocycle(edge,Ωij))
    end
    return cocycle
end

function sw1(cocycle::Vector{CechCocycle})
    #compute the first stiefel-whitney class of a given cechcocycle

    sw = Vector{Ripserer.ChainElement{Ripserer.Simplex{1,Float64,Int64},Mod{2}}}(undef,0)
    for elem in cocycle
        #add an 1 (mod 2) to the sw cochain if the transition function is -1 ((-1,1) maps to (0,1) as a Z/2 isomorphism)
        if det(transitionfunction(elem)) < 0
            chainelement = Ripserer.ChainElement{Ripserer.Simplex{1,Float64,Int64},Mod{2}}(cechsimplex(elem),1)
            push!(sw,chainelement)
        end
    end
    return Ripserer.Chain{Mod{2},Ripserer.Simplex{1,Float64,Int64}}(sw)
end
