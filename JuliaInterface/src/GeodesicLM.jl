module GeodesicLM

const geodesiclmlib = string(Pkg.dir(), "/GeodesicLM/builds/libgeodesiclm.so")


## Dictionary of convergence messages
converged_msg = { 1=>"artol", 2=>"CGoal", 3=>"gtol", 4=>"xtol", 5=>"xrtol", 6=>"ftol", 7=>"frtol", -1=>"iters", -2=>"nfev", -3=>"njev", -4=>"naev", -5=>"maxlam", -6=>"minlam", -10=>"user_termination", -11=>"func_fail"}

## Dummy functions
function jacobian_dummy(m_::Ptr{Int32}, n_::Ptr{Int32}, x_::Ptr{Float64}, fjac_::Ptr{Float64})
    return convert(Cint, 1)
end
    
function Avv_dummy(m_::Ptr{Int32}, n_::Ptr{Int32}, x_::Ptr{Float64}, v_::Ptr{Float64}, acc_::Ptr{Float64})
    return convert(Cint, 1)
end

function callback_dummy(m_::Ptr{Int32}, n_::Ptr{Int32}, x_::Ptr{Float64}, v_::Ptr{Float64}, a_::Ptr{Float64}, fvec_::Ptr{Float64}, fjac_::Ptr{Float64}, acc_::Ptr{Float64}, lam_::Ptr{Float64}, dtd_::Ptr{Float64}, fvec_new_::Ptr{Float64}, accepted_::Ptr{Int32}, info_::Ptr{Int32})
    return convert(Cint, 1)
end    

## Wrapper around the main geodesiclm fortran routine
function geodesiclm(func::Function, x::Array{Float64}, m::Integer, n::Integer;
                    jacobian::Function = jacobian_dummy, Avv::Function = Avv_dummy, callback::Function = callback_dummy,
                    center_diff::Bool = false, h1::Float64 = sqrt(eps()), h2::Float64 = eps()^0.25,
                    damp_mode::Integer = 1,
                    maxiter::Integer = 500, maxfev::Integer = 0, maxjev::Integer = 0, maxaev::Integer = 0, maxlam::Float64 = -1.0, minlam::Float64 = -1.0,
                    artol::Float64 = 1e-3, Cgoal::Float64 = sqrt(eps()), gtol::Float64 = sqrt(eps()), 
                    xtol::Float64 = sqrt(eps()), xrtol::Float64 = -1.0, ftol::Float64 = sqrt(eps()), frtol::Float64 = -1.0, 
                    print_level::Integer = 1, print_unit::Integer = 6,
                    imethod::Integer = 0, iaccel::Integer =1, ibold::Integer = 0, ibroyden::Integer = 0,
                    initialfactor::Float64 = 0.001, factoraccept::Float64 = 3.0, factorreject::Float64 = 2.0, avmax::Float64 = 0.75
                    )

    const func_f = cfunction(func, Cint, (Ptr{Int32}, Ptr{Int32}, Ptr{Float64}, Ptr{Float64}))
    const jacobian_f = cfunction(jacobian, Cint, (Ptr{Int32}, Ptr{Int32}, Ptr{Float64}, Ptr{Float64}))
    const Avv_f = cfunction(Avv, Cint, (Ptr{Int32}, Ptr{Int32}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}))
    const callback_f = cfunction(callback, Cint, (Ptr{Int32}, Ptr{Int32}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64},Ptr{Float64}, Ptr{Float64}, Ptr{Int32}, Ptr{Int32}) )
   
    
    analytic_jac = (jacobian != jacobian_dummy)
    analytic_avv = (Avv != Avv_dummy)

    fvec = zeros(m)
    fjac = zeros(m,n)

    # we need to get the value of info back from geodesiclm.  This seems to be the easiest way to make the changed value accessible after the call
    info = [convert(Int32,0)]


    dtd = eye(n)
    niters = [convert(Int32, 0)]
    nfev = [convert(Int32,0)]
    njev = [convert(Int32,0)]
    naev = [convert(Int32, 0)]
    converged = [convert(Int32,0)]
    
    ccall( (:geodesiclm_, geodesiclmlib), Void, 
          (Ptr{Void}, Ptr{Void}, Ptr{Void}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Int32}, Ptr{Int32}, Ptr{Void}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}), 
    func_f, jacobian_f, Avv_f, x, fvec, fjac, &convert(Int32, n), &convert(Int32, m), callback_f, info,
    &analytic_jac, &analytic_avv, &center_diff, &h1, &h2, dtd, &convert(Int32, damp_mode), 
    niters, nfev, njev, naev, 
    &convert(Int32, maxiter), &convert(Int32,maxfev), &convert(Int32,maxjev), &convert(Int32,maxaev),
    &maxlam, &minlam, &artol, &Cgoal, &gtol, &xtol, &xrtol, &ftol, &frtol, converged,
    &convert(Int32,print_level), &convert(Int32,print_unit),
    &convert(Int32,imethod), &convert(Int32,iaccel), &convert(Int32,ibold), &convert(Int32,ibroyden),
    &initialfactor, &factoraccept, &factorreject, &avmax)

    return x, {"converged"=>converged[1], "iters"=>[niters[1], nfev[1], njev[1], naev[1]], "msg"=>converged_msg[converged[1]], "fvec"=>fvec, "fjac"=>fjac}
end

export geodesiclm

end
