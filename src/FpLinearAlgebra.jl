using Ripserer
using Nemo


function chaintovector(complex::Ripserer.AbstractFiltration{<:Integer,<:Number},
                       chain::Ripserer.Chain{<:Integer,Ripserer.Simplex{1,Float64,Int64}},
                       p::Number)
        #send a (Ripserer) 1 chain to the associated vector modulo p which can be interpreted by Nemo
        #preferably p should be the same as that used to generate the chain (not strictly required when we want mod p reduction of a Z chain)

        field = GF(p)
        oneskeleton = complex[1]
        outputvector = zero_matrix(field,length(oneskeleton),1)

        #generate a lookup table for the simplex order
        oneskeletonlookup = Dict([(oneskeleton[i],i) for i in eachindex(oneskeleton)])
        for c in chain
            i = oneskeletonlookup[simplex(c)]
            coeff = Int(coefficient(c))
            outputvector[i] = field(coeff)
        end

        return outputvector
end