module GeodesicLM

## Dictionary of convergence messages
converged_msg = { 1=>"artol", 2=>"CGoal", 3=>"gtol", 4=>"xtol", 5=>"xrtol", 6=>"ftol", 7=>"frtol", -1=>"iters", -2=>"nfev", -3=>"njev", -4=>"naev", -5=>"maxlam", -6=>"minlam", -10=>"user_termination", -11=>"func_fail"}

## Dummy functions
function jacobian_dummy(m_::Ptr{Int64}, n_::Ptr{Int64}, x_::Ptr{Float64}, fjac_::Ptr{Float64})
    return convert(Cint, 1)
end
    
function Avv_dummy(m_::Ptr{Int64}, n_::Ptr{Int64}, x_::Ptr{Float64}, v_::Ptr{Float64}, acc_::Ptr{Float64})
    return convert(Cint, 1)
end

function callback_dummy(m_::Ptr{Int64}, n_::Ptr{Int64}, x_::Ptr{Float64}, v_::Ptr{Float64}, a_::Ptr{Float64}, fvec_::Ptr{Float64}, fjac_::Ptr{Float64}, acc_::Ptr{Float64}, lam_::Ptr{Float64}, dtd_::Ptr{Float64}, fvec_new_::Ptr{Float64}, accepted_::Ptr{Int64}, info_::Ptr{Int64})
    return convert(Cint, 1)
end    

## Wrapper around the main geodesiclm fortran routine
function geodesiclm(func::Function, x::Array{Float64}, m::Integer, n::Integer;
                    jacobian::Function = jacobian_dummy, Avv::Function = Avv_dummy, callback::Function = callback_dummy,
                    center_diff::Bool = false, h1::Float64 = sqrt(eps()), h2::Float64 = eps()^0.25,
                    damp_mode::Int64 = 1,
                    maxiter::Int64 = 500, maxfev::Int64 = 0, maxjev::Int64 = 0, maxaev::Int64 = 0, maxlam::Float64 = -1.0, minlam::Float64 = -1.0,
                    artol::Float64 = 1e-3, Cgoal::Float64 = sqrt(eps()), gtol::Float64 = sqrt(eps()), 
                    xtol::Float64 = sqrt(eps()), xrtol::Float64 = -1.0, ftol::Float64 = sqrt(eps()), frtol::Float64 = -1.0, 
                    print_level::Int64 = 1, print_unit::Int64 = 6,
                    imethod::Int64 = 0, iaccel::Int64 =1, ibold::Int64 = 0, ibroyden::Int64 = 0,
                    initialfactor::Float64 = 0.001, factoraccept::Float64 = 3.0, factorreject::Float64 = 2.0, avmax::Float64 = 0.75
                    )

    const func_f = cfunction(func, Cint, (Ptr{Int64}, Ptr{Int64}, Ptr{Float64}, Ptr{Float64}))
    const jacobian_f = cfunction(jacobian, Cint, (Ptr{Int64}, Ptr{Int64}, Ptr{Float64}, Ptr{Float64}))
    const Avv_f = cfunction(Avv, Cint, (Ptr{Int64}, Ptr{Int64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}))
    const callback_f = cfunction(callback, Cint, (Ptr{Int64}, Ptr{Int64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64},Ptr{Float64}, Ptr{Float64}, Ptr{Int64}, Ptr{Int64}) )
   

    analytic_jac = (jacobian != jacobian_dummy)
    analytic_avv = (Avv != Avv_dummy)

    fvec = zeros(m)
    fjac = zeros(m,n)

    # we need to get the value of info back from geodesiclm.  This seems to be the easiest way to make the changed value accessible after the call
    info = [0]


    dtd = eye(n)
    niters = [0]
    nfev = [0]
    njev = [0]
    naev = [0]
    converged = [0]
    
    ccall( (:geodesiclm_, "./libgeodesiclm.so"), Void, 
          (Ptr{Void}, Ptr{Void}, Ptr{Void}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Int64}, Ptr{Int64}, Ptr{Void}, Ptr{Int64}, Ptr{Int64}, Ptr{Int64}, Ptr{Int64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Int64}, Ptr{Int64}, Ptr{Int64}, Ptr{Int64}, Ptr{Int64}, Ptr{Int64}, Ptr{Int64}, Ptr{Int64}, Ptr{Int64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Int64}, Ptr{Int64}, Ptr{Int64}, Ptr{Int64}, Ptr{Int64}, Ptr{Int64}, Ptr{Int64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}), 
    func_f, jacobian_f, Avv_f, x, fvec, fjac, &n, &m, callback_f, info, &analytic_jac, &analytic_avv, &center_diff, &h1, &h2, dtd, &damp_mode, 
    niters, nfev, njev, naev, 
    &maxiter, &maxfev, &maxjev, &maxaev, &maxlam, &minlam, &artol, &Cgoal, &gtol, &xtol, &xrtol, &ftol, &frtol, converged, &print_level, &print_unit, &imethod, &iaccel, &ibold, &ibroyden, &initialfactor, &factoraccept, &factorreject, &avmax)

    return x, {"converged"=>converged[1], "iters"=>[niters[1], nfev[1], njev[1], naev[1]], "msg"=>converged_msg[converged[1]], "fvec"=>fvec, "fjac"=>fjac}
end

export geodesiclm

end
