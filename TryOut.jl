include("MGO.jl")
using ForwardDiff
using OrdinaryDiffEq
using Plots

function D(x::AbstractVector,k::AbstractVector)
    -k.*k - x
end

τmin = 0
τsteps = 1000
τmax = 1
τ = range(τmin,τmax,1000)
sol = TraceRayIP([1.0,1.0,0,0,0,1.1],0,τmin,τmax,D::Function)
plot(sol,show=true)
