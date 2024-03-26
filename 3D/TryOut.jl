include("MGO.jl")
using ForwardDiff
using OrdinaryDiffEq
using Plots

#Depends on system
function D(x::AbstractVector,k::AbstractVector)
    -k.*k - x
end

function ComputeRayManifold(τmin,τmax,kx0,ky0,D::Function,y0s::AbstractVector)
    xk0 = [0.0,y0s[1],0,kx0,ky0,0.1]
    t,x,y,z,kx,ky,kz = TraceRayIP(xk0,0,τmin,τmax,D::Function)
    Ray = [t,x,y,z,kx,ky,kz]
    RayManifold = [Ray]
    for y0 in y0s[2:size(y0s)[1]-1]
        xk0 = [0.0,y0,0,kx0,ky0,0.1]
        t,x,y,z,kx,ky,kz = TraceRayIP(xk0,0,τmin,τmax,D::Function)
        Ray = [t,x,y,z,kx,ky,kz]
        RayManifold = [RayManifold;[Ray[:]]]
    end
    RayManifold
end

τmin = 0
τsteps = 1000
τmax = 1
τ = range(τmin,τmax,1000)
#xk0 = [1.0,1.0,0,0,0,1.1]
#t,x,y,z,kx,ky,kz = TraceRayIP(xk0,0,τmin,τmax,D::Function)
y0s = LinRange(-1,1,20)
RayManifold = ComputeRayManifold(τmin,τmax,0.01,0.7,D::Function,y0s)
plot()
for Ray in RayManifold
    plot!(Ray[2],Ray[3],show=true)
end
