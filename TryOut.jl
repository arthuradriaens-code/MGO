include("MGO.jl")
using ForwardDiff
using OrdinaryDiffEq
using Plots

function D(x::AbstractVector,k::AbstractVector)
    -k.*k - x
end

#τmin = zeros(1,6)
τmin = 0
#τmax = ones(1,6)
τmax = 1
sol = TraceRay([1.0,1.0,0,0,0,1.1],0,τmin,τmax,D::Function)
plot(sol)
