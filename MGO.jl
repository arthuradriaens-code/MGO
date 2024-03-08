using ForwardDiff
using OrdinaryDiffEq

function TraceRay(xk0::AbstractVector,ω0,τmin,τmax,D,rtol=1e-3,τsteps=1000)
    # xk = [x,y,z,kx,ky,kz]
    function RHS(xk::AbstractVector,p,t)
        x_ = xk[1:3]
        k_ = xk[4:6]
        ∂D∂k(x,k) = ForwardDiff.gradient(D(x,k),k)
        ∂D∂x(x,k) = ForwardDiff.gradient(D(x,k),x)
        print(∂D∂x(x,k))
        [-∂D∂k,∂D∂x]
    end

    prob = ODEProblem(RHS,xk0,(τmin,τmax))
    sol = solve(prob,Tsit5(),reltol=rtol)
end
