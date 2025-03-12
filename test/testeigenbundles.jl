using Random
using LinearAlgebra
using Ripserer
using DiscreteVectorBundles

function sphere(n)
    values = rand(Xoshiro(0),Float64,(n,2))
    values[:,2] = 2*π*values[:,2]
    values[:,1] = 2*values[:,1]-ones(n)
    values[:,1] = map(acos,values[:,1])

    pointcloud = [(sin(values[i,1])*cos(values[i,2]),sin(values[i,1])*sin(values[i,2]),cos(values[i,1])) for i in 1:n]
    return pointcloud
end

function annulus(n, r1=1, r2=1.1, offset=(0, 0))
    result = Tuple{Float64,Float64}[]
    while length(result) < n
        point = 2 * r2 * rand(2) .- r2
        if r1 < norm(point) < r2
            push!(result, (point[1] + offset[1], point[2] + offset[2]))
        end
    end
    return result
end

function H_mobius(n)
    return Hermitian([n[1] n[2];
                     n[2] -n[1]]) #making sure the matrix is strictly hermitian
end

function H_sphere(n)
    nvec = [n[1] n[2] n[3]]
    return Hermitian(I(3)-nvec'*nvec)
end

function H_geofluid(n)
    return Hermitian([
        [0 n[1] n[2]];
        [n[1] 0 n[3]*1im];
        [n[2] -n[3]*1im 0]
    ])
end

function mobiusbundletest()
    N=200
    mobius_data = annulus(N)
    localtriv = constructrealeigenbundle(mobius_data,H_mobius,[2])
    dth = approxcocycledeath(localtriv)-1e-2
    a,cocyclebasis = solveincohomologybasis(Alpha(mobius_data),sw1(localtriv,dth),2,dth)
    return a[1]
end

function spherebundletest()
    #chech that the euler class of the tangent bundle to the sphere is calculated correctly
    N=200
    spherepoints = sphere(N)

    spherebundle = [eigvecs(H_sphere(spherepoints[i])) for i in 1:N]
    tangentbundle = [spherebundle[i][:,2:3] for i in 1:N]
    sphere_alpha_complex = Alpha(spherepoints)

    local_triv = LocalTriv(sphere_alpha_complex,tangentbundle)
    dth = approxcocycledeath(local_triv,1)-1e-2
    orientbundle!(local_triv,dth)

    v1 = local_triv.bases[1][:,1]
    v2 = local_triv.bases[1][:,2]

    if dot(cross(v1,v2),spherepoints[1]) < 0
        fliporientation!(local_triv)
    end

    p = 3
    c1 = chaintovector(sphere_alpha_complex,eu(local_triv,dth),p)
    μS2 = chaintovector(sphere_alpha_complex,μ(spherepoints,dth,p),p)
    return evaluate(μS2,c1)
end

function geofluidstest(branch)
    #chec that the chern numbers of the Weyl point in the f plane model (from geophysical fluid dynamics) are computed correctly

    N = 200
    spherepoints = sphere(N)

    localtriv = constructcomplexeigenbundle(spherepoints,H_geofluid,branch)
    dth = approxcocycledeath(localtriv,1)-1e-2
    orientbundle!(localtriv,dth)

    p = 3
    c1 = chaintovector(filtration(localtriv),eu(localtriv,dth),p)
    μS2 = chaintovector(filtration(localtriv),μ(spherepoints,dth,p),p)
    return evaluate(μS2,c1)
end