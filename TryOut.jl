include("MGO.jl")
using ForwardDiff
using OrdinaryDiffEq
using Plots

#Depends on system
function D(x::AbstractVector,k::AbstractVector)
    -k.*k - x
end

τmin = 0
τsteps = 1000
τmax = 1
τ = range(τmin,τmax,1000)
xk0 = [1.0,1.0,0,0,0,1.1]
t,x,y,z,kx,ky,kz = TraceRayIP(xk0,0,τmin,τmax,D::Function)
plot(x,y,z,show=true)
