using Ripserer

function geometricorientation(pointcloud::Vector{Tuple{Float64, Float64, Float64}},
                              complex::Ripserer.AbstractFiltration{<:Integer,<:Number},
                              testtriangle::Ripserer.Simplex{2,Float64,Int64},
                              threshold::Number)
        #Take a pointcloud, a complex, and a triangle in the complex and determine a chainelement 
        #giving the geometric orientation (outwards pointing) induced by the point cloud embedding assuming
        #the complex is topologically a connected compact 2-manifold
        
        #Works by shooting a ray in the direction of the normal (from the right hand rule) to a single triangle
        #and counting how many times it intersects to complex. If this intersection number is odd the triangle is 
        #inwards pointing, if it is even it is an outwards pointing triangle.

        v1,v2,v3 = vertices(testtriangle)
        n = cross(collect(pointcloud[v2])-collect(pointcloud[v1]),
                  collect(pointcloud[v3])-collect(pointcloud[v2]))
        b = collect(pointcloud[v2])
        bindex = v2 

        intersectionnumber = 0
        #shoot a half ray t*n+b (t>0) and count the number of intersections with triangles in the complex (at death time = threshold)
        for triangle in complex[2]
            if triangle != testtriangle && birth(triangle) <= threshold
                #each triangle specifies a plane by its normal and a basepoint. P = {x: mᵀ(x-basepoint) =0}
                v1,v2,v3 = vertices(triangle)

                if any([v1==bindex,v2==bindex,v3==bindex]) #ignore triangles which share the vertex the ray is based at (false intersections)
                    continue
                end

                basepoint = collect(pointcloud[v2])
                e1 = collect(pointcloud[v1])-collect(pointcloud[v2]) #basis for the plane in which the triangle lines
                e2 = collect(pointcloud[v3])-collect(pointcloud[v2])
                m = cross(-e1,e2)

                mdotn = dot(m,n)
                if mdotn == 0
                    continue
                else
                    t = dot(m,basepoint-b)/mdotn #solve for the intersection point on the ray
                    if t > 0 #check if the intersection is on the positive side of the ray
                        y = t*n+b
                
                        M = [e1 e2]
                        α,β  = M\(y-basepoint) #solve for the intersection point in the basis of the triangle edges

                        if α >= 0 && β >= 0 && α+β <=1
                            intersectionnumber += 1
                        end
                    end
                end
            end
        end
        if mod(intersectionnumber,2) == 0
            return 1
        else
            return -1
        end
    end

function μ(pointcloud::Vector{Tuple{Float64, Float64, Float64}},
           threshold::Number,
           p::Number)
        #pull the fundamental homology class of the Alpha filtration of a point cloud mod p
        #with the orientation chosen so it agrees with the outward pointing vector field
        #assuming that poincloud is sampled from a connected compact surface in R³

        complex = Alpha(pointcloud)
        diagram = ripserer(Alpha,pointcloud,modulus=p,dim_max = 2; alg= :involuted)
        fundamentalclass = representative(diagram[3][end])

        
        #take a simplex in the topological fundamental class and compute its geometric orientation
        testtriangle = simplex(fundamentalclass[1])
        if Mod{p}(geometricorientation(pointcloud,complex,testtriangle,threshold)) == coefficient(fundamentalclass[1])
            return fundamentalclass
        elseif Mod{p}(geometricorientation(pointcloud,complex,testtriangle,threshold)) == -coefficient(fundamentalclass[1])
            fundamentalclass2 = [Ripserer.ChainElement{Ripserer.Simplex{2,Float64,Int64},Mod{p}}(simplex(triangle),-coefficient(triangle)) for triangle in fundamentalclass]
            return Ripserer.Chain{Mod{p},Ripserer.Simplex{2,Float64,Int64}}(fundamentalclass2)
        else
            throw("H2 chain has elements which are not 1 or -1, point cloud may be faulty") 
        end
end

function constructrealeigenbundle(pointcloud::Vector{<:Tuple},
                              H::Function,
                              branches::Vector{Int})
        #Take a pointcloud and a self-adjoint operator and construct a local trivilisation of the discrete vector bundle 
        #over the alpha filtration of the pointcloud. 

        #Note that the eigenvectors of H will be calculated in the order indexed by their eigenvalue. 
        #The ``branches'' input is a tuple {For example/ [1] or [2,3]} which specifies which eigenvectors 
        #to include int the bundle by their eigenvalue order.

        N = length(pointcloud)
        complex = Alpha(pointcloud)
        trivialbundle = [eigvecs(H(pointcloud[i])) for i in 1:N]
        eigenbundle = [trivialbundle[i][:,branches] for i in 1:N]

        return LocalTriv(complex,eigenbundle)
end

function constructcomplexeigenbundle(pointcloud::Vector{<:Tuple},
                                     H::Function,
                                     branch::Int)
        #Take a pointcloud and a self-adjoint operator and construct a local trivilisation of a discrete complex line bundle 
        #over the alpha filtration of the pointcloud returned as a real line bundle with the canonical orientation
        #induced by the complex structure.

        #Note that the eigenvectors of H will be calculated in the order indexed by their eigenvalue. 
        #The ``branch'' input is a integer which specifies which eigenvalue we want the eigenvectors for.

        N = length(pointcloud)
        complex = Alpha(pointcloud)
        trivialbundle = [eigvecs(H(pointcloud[i])) for i in 1:N]
        eigenbundle = [canonicallyorientatedplane(trivialbundle[i][:,branch]) for i in 1:N]

        return LocalTriv(complex,eigenbundle)
end

#Core functions describing the isomorphism from complex vector bundles to oriented real bundles
function j(n)
    #2nx2n matrix specifying the almost complex structure
    return [zeros(n,n) -I(n); I(n) zeros(n,n)]
end

function comptoreal(complexvector)
    #map a complex n vector into R^n
    return [real(complexvector);imag(complexvector)]
end

function canonicallyorientatedplane(complexvector)
    #take a complex vector and construct the real plane it defines (as a Stiefel matrix)
    #given the canonical orientation induced by multiplication by j

    realvec = comptoreal(complexvector)
    return [realvec j(length(complexvector))*realvec]
end