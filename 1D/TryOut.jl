include("MGO.jl")
using ForwardDiff
using OrdinaryDiffEq
using SpecialFunctions
using Plots

#Depends on system
function D(x,k)
    -x - k*k
end

τ1min = -8
τ1steps = 1000
τ1max = 0
τ1 = range(τ1min,τ1max,τ1steps)
plot(τ1,airyai.(τ1),show=true) #reference exact solution
ϕ0 = Ai(τ1min)
