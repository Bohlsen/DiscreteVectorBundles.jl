using Ripserer
using Nemo


function chaintovector(complex::Ripserer.AbstractFiltration{<:Integer,<:Number},
                       chain::Ripserer.Chain{<:Integer,Ripserer.Simplex{1,Float64,Int64}},
                       p::Number,
                       threshold::Number = Inf)
        #send a (Ripserer) 1 chain to the associated vector modulo p which can be interpreted by Nemo
        #preferably p should be the same as that used to generate the chain (not strictly required when we want mod p reduction of a Z chain)

        field, = residue_ring(ZZ,p)

        oneskeleton = complex[1]
        outputvector = matrix(field,zero_matrix(ZZ,length(oneskeleton),1))

        #generate a lookup table for the simplex order
        oneskeletonlookup = Dict([(oneskeleton[i],i) for i in eachindex(oneskeleton)])
        for c in chain
            if birth(simplex(c)) > threshold
                continue
            end
            i = oneskeletonlookup[simplex(c)]
            coeff = Int(coefficient(c))
            outputvector[i] = field(coeff)
        end

        return outputvector
end

function chaintovector(complex::Ripserer.AbstractFiltration{<:Integer,<:Number},
                       chain::Ripserer.Chain{<:Integer,Ripserer.Simplex{2,Float64,Int64}},
                       p::Number,
                       threshold::Number = Inf)
                        #send a (Ripserer) 2 chain to the associated vector modulo p which can be interpreted by Nemo
                        #preferably p should be the same as that used to generate the chain (not strictly required when we want mod p reduction of a Z chain)

        field, = residue_ring(ZZ,p)

        twoskeleton = complex[2]
        outputvector = matrix(field,zero_matrix(ZZ,length(twoskeleton),1))

        #generate a lookup table for the simplex order
        twoskeletonlookup = Dict([(twoskeleton[i],i) for i in eachindex(twoskeleton)])
        for c in chain
            if birth(simplex(c)) >threshold
                continue
            end
            i = twoskeletonlookup[simplex(c)]
            coeff = Int(coefficient(c))
            outputvector[i] = field(coeff)
        end

        return outputvector
end

function one∂(complex::Ripserer.AbstractFiltration{<:Integer,<:Number},
              p::Number,
              threshold::Number=Inf)
        #create the ∂ matrix for 1 chains from a filtration only including edges with birth time less than threshold

        field, = residue_ring(ZZ,p)

        oneskeleton = complex[1]
        zeroskeleton = complex[0]

        outputmatrix = matrix(field,zero_matrix(ZZ,length(zeroskeleton),0))

        #generate a lookup table for the simplex order
        zeroskeletonlookup = Dict([(zeroskeleton[i],i) for i in eachindex(zeroskeleton)])

        for edge in oneskeleton
            if birth(edge) <= threshold
                v1,v2 = Ripserer.boundary(complex,edge) #pull the boundary of the edge
                gen = matrix(field,zero_matrix(ZZ,length(zeroskeleton),1))

                gen[zeroskeletonlookup[v2]] = -1
                gen[zeroskeletonlookup[v1]] = 1

                outputmatrix = matrix(field,[outputmatrix gen])
            else
                #include the extra edges but as zero rows so they pad the vector but don't add to the column space
                gen = zero_matrix(field,length(zeroskeleton),1)
                outputmatrix = matrix(field,[outputmatrix gen])
            end
        end
        return outputmatrix
end

function two∂(complex::Ripserer.AbstractFiltration{<:Integer,<:Number},
              p::Number,
              threshold::Number=Inf)
        field, = residue_ring(ZZ,p)

        oneskeleton = complex[1]
        twoskeleton = complex[2]

        outputmatrix = matrix(field,zero_matrix(ZZ,length(oneskeleton),0))

        #generate a lookup table for the simplex order
        oneskeletonlookup = Dict([(oneskeleton[i],i) for i in eachindex(oneskeleton)])

        for triangle in twoskeleton
            if birth(triangle) <= threshold
                v1,v2,v3 = Ripserer.boundary(complex,triangle) #pull the boundary of the edge
                gen = matrix(field,zero_matrix(ZZ,length(oneskeleton),1))

                gen[oneskeletonlookup[v3]] = 1
                gen[oneskeletonlookup[v2]] = -1
                gen[oneskeletonlookup[v1]] = 1

                outputmatrix = matrix(field,[outputmatrix gen])
            else
                #include the extra edges but as zero rows so they pad the vector but don't add to the column space
                gen = zero_matrix(field,length(oneskeleton),1)
                outputmatrix = matrix(field,[outputmatrix gen])
            end
        end
        return outputmatrix
end


function zeroδ(complex::Ripserer.AbstractFiltration{<:Integer,<:Number},
                p::Number,
                threshold::Number=Inf)
        #compute the zero coboundary matrix as the transpose of the one∂
        return transpose(one∂(complex,p,threshold))
end

function oneδ(complex::Ripserer.AbstractFiltration{<:Integer,<:Number},
                p::Number,
                threshold::Number=Inf)
        #compute the two coboundary matrix as the transpose of the two∂
        return transpose(two∂(complex,p,threshold))
end


function findpivot(row::Vector{zzModRingElem})
    for i in eachindex(row)
        if row[i] != 0
            return i
        end
    end
    return 0
end


function linearsolvemodp(M::zzModMatrix,b::zzModMatrix,p::Number)
    #Solve M*a=b over the Fp field

    R, = residue_ring(ZZ,p)
    A = matrix(R,[M b]) #construct the augmented matrix
    Aref = strong_echelon_form(A) #compute the row echelon form of the augmented matrix


    Mref = Aref[:,1:end-1]
    bref = Aref[:,end]

    a = matrix(R,zero_matrix(ZZ,size(Mref)[2],1))
    pivots = [findpivot(Mref[i,:]) for i in 1:size(Mref)[1]]

    for i in reverse(eachindex(pivots)) #perform backsubstitution to solve for the a vector
        if pivots[i] == 0
            if bref[i] != 0 
                throw("nonzero b in zero row, cannot find a solution")
            else
                continue
            end
        else
            a[pivots[i]] = bref[i]
            for j in 1:i-1
                bref[j] = bref[j] - Mref[j,i]*bref[i]
            end
        end
    end
    return a
end

function evaluate(chain::zzModMatrix,cochain::zzModMatrix)
    #evaluate a cochain against a specific chain
    if length(chain) != length(cochain)
        throw("Lengths must be equal")
    elseif modulus(chain[1]) != modulus(cochain[1])
        throw("chain and cochain must be over the same ring")
    end
    return (transpose(cochain)*chain)[1]
end