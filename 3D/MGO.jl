using ForwardDiff
using OrdinaryDiffEq

function TraceRayIP(xk0,ω0,τmin,τmax,D,τsteps=1000)
    # xk = [x,y,z,kx,ky,kz]
    function RHS!(dxk,xk,p,t) #Ray Hamilton's equations
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
    probf = ODEProblem(RHS!,xk0,(τmin,τmax)) #forward in time
    solf = solve(probf,Tsit5(),saveat=0.1)
    probb = ODEProblem(RHS!,xk0,(τmin,τmin-0.2)) #back in time
    solb = solve(probb,Tsit5(),saveat=0.1)

    #joining backwards and forwards, defining t,x,y,z,kx,ky,kz
    x = [reverse(reduce(vcat,transpose.(solb.u))[:,1]);reduce(vcat,transpose.(solf.u))[:,1]]
    y = [reverse(reduce(vcat,transpose.(solb.u))[:,2]);reduce(vcat,transpose.(solf.u))[:,2]]
    z = [reverse(reduce(vcat,transpose.(solb.u))[:,3]);reduce(vcat,transpose.(solf.u))[:,3]]
    kx = [reverse(reduce(vcat,transpose.(solb.u))[:,4]);reduce(vcat,transpose.(solf.u))[:,4]]
    ky = [reverse(reduce(vcat,transpose.(solb.u))[:,5]);reduce(vcat,transpose.(solf.u))[:,5]]
    kz = [reverse(reduce(vcat,transpose.(solb.u))[:,6]);reduce(vcat,transpose.(solf.u))[:,6]]
    t = [reverse(reduce(vcat,transpose.(solb.t)));reduce(vcat,transpose.(solf.t))]

    return t,x,y,z,kx,ky,kz
end
function TraceRayIP2D(xk0,ω0,τmin,τmax,D,τsteps=1000)
    # xk = [x,y,kx,ky]
    function RHS!(dxk,xk,p,t) #Ray Hamilton's equations
        x_ = xk[1:2]
        k_ = xk[3:4]
        JD(xk) = ForwardDiff.jacobian(xk -> D(x_,k_),xk)
        dxk[1] = JD(xk)[1,4]
        dxk[2] = JD(xk)[2,4]
        dxk[3] = -JD(xk)[1,1]
        dxk[4] = -JD(xk)[2,2]
    end
    probf = ODEProblem(RHS!,xk0,(τmin,τmax)) #forward in time
    solf = solve(probf,Tsit5(),saveat=0.1)
    probb = ODEProblem(RHS!,xk0,(τmin,τmin-0.2)) #back in time
    solb = solve(probb,Tsit5(),saveat=0.1)

    #joining backwards and forwards, defining t,x,y,z,kx,ky,kz
    x = [reverse(reduce(vcat,transpose.(solb.u))[:,1]);reduce(vcat,transpose.(solf.u))[:,1]]
    y = [reverse(reduce(vcat,transpose.(solb.u))[:,2]);reduce(vcat,transpose.(solf.u))[:,2]]
    kx = [reverse(reduce(vcat,transpose.(solb.u))[:,3]);reduce(vcat,transpose.(solf.u))[:,3]]
    ky = [reverse(reduce(vcat,transpose.(solb.u))[:,4]);reduce(vcat,transpose.(solf.u))[:,4]]
    t = [reverse(reduce(vcat,transpose.(solb.t)));reduce(vcat,transpose.(solf.t))]
    return t,x,y,kx,ky
end
