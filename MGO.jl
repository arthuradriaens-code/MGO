using ForwardDiff
using OrdinaryDiffEq

function TraceRayIP(xk0,ω0,τmin,τmax,D,τsteps=1000)
    # xk = [x,y,z,kx,ky,kz]
    function RHS!(dxk,xk,p,t) #Ray Hamilton's equationS
        x_ = xk[1:3]
        k_ = xk[4:6]
        JD(xk) = ForwardDiff.jacobian(xk -> D(xk[1:3],xk[4:6]),xk)
        dxk[1] = JD(xk)[1,4]
        dxk[2] = JD(xk)[2,5]
        dxk[3] = JD(xk)[3,6]
        dxk[4] = -JD(xk)[1,1]
        dxk[5] = -JD(xk)[2,2]
        dxk[6] = -JD(xk)[3,3]
    end
    prob = ODEProblem(RHS!,xk0,(τmin,τmax))
    sol = solve(prob,Tsit5())
end
