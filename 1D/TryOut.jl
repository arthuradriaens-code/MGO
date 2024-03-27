include("MGO.jl")
using ForwardDiff
using OrdinaryDiffEq
using SpecialFunctions
using Plots

#Depends on system
function D(xk::AbstractVector)
    x = xk[1]
    k = xk[2]
    -x - k*k
end

τ1min = 0
τ1steps = 1000
τ1max = 8
τ1 = range(τ1min,τ1max,τ1steps)
#plot(τ1,airyai.(τ1),show=true) #reference exact solution
xk0 = [-8,sqrt(8)]
phi0 = airyai(-8)
sol,t,nt,xs,ks,zs =  TraceRay(xk0,τ1min,τ1max,D::Function,-8)
itpx,itpk = get_go_field(sol,phi0)
plot(xs,ks,show=true)
#ϕ0 = airyai(τ1min)
