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

function H_sphere(n)
    nvec = [n[1] n[2] n[3]]
    return Hermitian(I(3)-nvec'*nvec)
end

function spherebundletest()
    N=200
    spherepoints = sphere(N)
    X = [spherepoints[i][1] for i in 1:N]
    Y = [spherepoints[i][2] for i in 1:N]
    Z = [spherepoints[i][3] for i in 1:N]

    spherebundle = [eigvecs(H_sphere(spherepoints[i])) for i in 1:N]
    tangentbundle = [spherebundle[i][:,2:3] for i in 1:N]
    sphere_alpha_complex = Alpha(spherepoints)

    local_triv = LocalTriv(sphere_alpha_complex,tangentbundle)
    dth = approxcocycledeath(local_triv,1)-1e-15
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