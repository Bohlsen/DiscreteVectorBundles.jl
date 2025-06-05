using Random
using LinearAlgebra
using Ripserer
using DiscreteVectorBundles

#pointcloud generating functions
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

#hermitian operators to compute the characteristic classes of
function H_mobius(n)
    return Hermitian([n[1] n[2];
                     n[2] -n[1]]) #making sure the matrix is strictly hermitian
end

function H_sphere(n)
    nvec = [n[1] n[2] n[3]]
    return Hermitian(I(3)-nvec'*nvec)
end

function H_geofluid(n)
    f = n[1]
    kx = n[2]
    ky = n[3]
    return Hermitian([
        [0 -f*1im kx];
        [f*1im 0 ky];
        [kx ky 0]
    ])
end

function H_dirac(n)
    px = n[1]
    py = n[2]
    m = n[3]

    return Hermitian([
        [m  px-py*1im];
        [px+py*1im -m]])
end

function H_TLCW(n,kz)
    ωp = n[1]
    kx = n[2]
    ky = n[3]
   
    ezcross = [[0 -1 0];
               [1 0 0];
               [0 0 0]]
    kcross = [[0 -kz ky];
               [kz 0 -kx];
               [-ky kx 0]]

    Hr = [[0*I(3) 0*I(3) 0*I(3)];
          [0*I(3) 0*I(3) -kcross];
          [0*I(3) kcross 0*I(3)]]

    Hi = [[ezcross -ωp*I(3) 0*I(3)];
          [ωp*I(3) 0*I(3) 0*I(3)];
          [0*I(3) 0*I(3) 0*I(3)]]

    return Hr+Hi*1im
end

#Helper functions for testing the extraction of the monopole index from synthetic experimental mobius_data
function Laguerrel(n,x)
    l=0
    for k in 0:n
        l += ((-1)^k*binomial(n,k)/factorial(k))*x^k
    end
    return l
end


#tests of the discrete characterstic class calculation
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

function diractest()
        N = 200
        spherepoints = sphere(N)
        function orientedlocaltriv_dirac(branch)
            localtriv = constructcomplexeigenbundle(spherepoints,H_dirac,branch)
            dth = approxcocycledeath(localtriv,1)-1e-2
            orientbundle!(localtriv,dth)
            return localtriv,dth
        end

        bundles = [orientedlocaltriv_dirac(branch) for branch in 1:2]

        chernnumbers = []
        p = 3

        for bundle in bundles
            localtriv = bundle[1]
            dth = bundle[2]
            c1 = chaintovector(filtration(localtriv),eu(localtriv,dth),p)
            μS2 = chaintovector(filtration(localtriv),μ(spherepoints,dth,p),p)
            push!(chernnumbers,evaluate(μS2,c1))
        end

        return chernnumbers
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

function TLCWtest(branch)
    #THIS IS CURRENTLY CONSIDERED EXPERIMENTAL as the calculation is consistent
    #but disagrees with Fu & Qin. 
    
    #It is possible that Fu & Qin contains a sign error (or there exists)
    #a subtle in the construction of the index theorem (a detailed investigation here is needed)

    ωpc = 0.58
    f = kz -> (sqrt(kz^4+4*kz^2)-kz^2)/2-ωpc

    N = 50 #intentially low point count here
    spherepoints = sphere(N)
    spherepoints = [(0.1*point[1]+ωpc,0.1*point[2],0.1*point[3]) for point in spherepoints]

    function bisection(f::Function,x0i::Number,x1i::Number,ε::Number)
        #find a zero of a single variable real function by the bisection method.
        #Either runs for 10000 bisections or if returns the current value if it is less that ε

        #assumes the x0i and x1i are given on the negative and positive sides of the zero respectively

        #NOTE THIS IS JUST FOR ONE TEST AND SHOULD NOT BE USED GENERALLY

        x0 = x0i
        x1 = x1i

        for i in 1:10000
            x = (x0+x1)/2
            if abs(f(x))<ε
                return x
            elseif f(x) < 0
                x0 = x
            else
                x1 = x
            end
        end
        throw("Did not find root")
    end

    kz = bisection(f,0.5,1,1e-10)

    localtriv = constructcomplexeigenbundle(spherepoints,n->H_TLCW(n,kz),branch)
    dth = approxcocycledeath(localtriv,1)-1e-2
    orientbundle!(localtriv,dth)

    p = 5

    c1 = chaintovector(filtration(localtriv),eu(localtriv,dth),p)
    μS2 = chaintovector(filtration(localtriv),μ(spherepoints,dth,p),p)
    return evaluate(μS2,c1)
end

