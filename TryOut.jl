include("MGO.jl")
using ForwardDiff
using OrdinaryDiffEq

function D(x::AbstractVector,k::AbstractVector)
    -k.*k - x
end
function D(x::Float64,k::Float64)
    -k.*k - x
end

sol = TraceRay([1.0,1.0,0,0,0,1.1],0,0,10,D::Function)
print(sol)
