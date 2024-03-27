using ForwardDiff
using OrdinaryDiffEq
using Interpolations

function TraceRay(xk0,τmin,τmax,D,boundaryx,τsteps=1000)
    # terminate integration if boundary is reached
    condition(xk, t, integrator) = xk[1] - boundaryx
    affect!(integrator) = terminate!(integrator)
    cb = ContinuousCallback(condition, affect!)

    function RHS!(dxk,xk,p,t) #Ray Hamilton's equations
        JD(xk) = ForwardDiff.gradient(xk -> D(xk),xk)
        ∂xD = JD(xk)[1]
        ∂kD = JD(xk)[2]
        dxk[1] = -∂kD
        dxk[2] = ∂xD
    end
    prob = ODEProblem(RHS!,xk0,(τmin,τmax)) #forward in time
    sol = solve(prob,Tsit5(),saveat=0.1,callback=cb)
    #probb = ODEProblem(RHS!,xk0,(τmin,τmin-0.2)) #back in time
    #solb = solve(probb,Tsit5(),saveat=0.1)

    #joining backwards and forwards, defining t,x,y,z,kx,ky,kz
    #x = [reverse(reduce(vcat,transpose.(solb.u))[:,1]);reduce(vcat,transpose.(solf.u))[:,1]]
    #kx = [reverse(reduce(vcat,transpose.(solb.u))[:,4]);reduce(vcat,transpose.(solf.u))[:,4]]
    #t = [reverse(reduce(vcat,transpose.(solb.t)));reduce(vcat,transpose.(solf.t))]
    extract_sol(sol)
end
function extract_sol(sol)
    xs = first.(sol.u);
    ks = last.(sol.u);
    #create matrix
    zs = hcat(xs, ks);
    t = sol.t;
    nt = size(t);
    return sol,t,nt,xs,ks,zs
end
function get_go_field(sol,phi0)
    sol,t,nt,xs,ks,zs = extract_sol(sol);
    #https://discourse.julialang.org/t/
    #differentiation-without-explicit-function-np-gradient/57784/2
    itpx = interpolate((t,), xs, Gridded(Linear()));
    itpk = interpolate((t,), ks, Gridded(Linear()));
    ∂τx = Interpolations.gradient.(Ref(itpx),t);
    ∂τk = Interpolations.gradient.(Ref(itpk),t);
    return itpx,itpk
end
