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
solf = TraceRayIP(xk0,0,τmin,τmax,D::Function)
solb = TraceRayIP(xk0,0,τmin,τmin-0.2,D::Function)
plot(solf,show=true)
