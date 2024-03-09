using ForwardDiff
using OrdinaryDiffEq

function TraceRay(xk0::AbstractVector,ω0,τmin,τmax,D,rtol=1e-3,τsteps=1000)
    # xk = [x,y,z,kx,ky,kz]
    function RHS(xk::AbstractVector,p,t)
        x_ = xk[1:3]
        k_ = xk[4:6]
        JD(xk) = ForwardDiff.jacobian(xk -> D(xk[1:3],xk[4:6]),xk)
        [JD(xk)[1,1],JD(xk)[2,2],JD(xk)[3,3],JD(xk)[1,4],JD(xk)[2,5],JD(xk)[3,6]]
    end

    prob = ODEProblem(RHS,xk0,(τmin,τmax))
    sol = solve(prob,Tsit5(),reltol=rtol)
end
