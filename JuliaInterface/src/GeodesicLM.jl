module GeodesicLM

const geodesiclmlib = string(Pkg.dir(), "/GeodesicLM/builds/libgeodesiclmwrapper.so")

## Dictionary of convergence messages
converged_msg = { 1=>"artol", 2=>"CGoal", 3=>"gtol", 4=>"xtol", 5=>"xrtol", 6=>"ftol", 7=>"frtol", -1=>"iters", -2=>"nfev", -3=>"njev", -4=>"naev", -5=>"maxlam", -6=>"minlam", -10=>"user_termination", -11=>"func_fail"}



## c functions
function func_(m_::Ptr{Int32}, n_::Ptr{Int32}, x_::Ptr{Float64}, fvec_::Ptr{Float64}, f_::Ptr{Void})
    m = unsafe_load(m_)
    n = unsafe_load(n_)
    x = pointer_to_array(x_, (n,))
    fvec = pointer_to_array(fvec_, (m,))
    f = unsafe_pointer_to_objref(f_)::Function
    fvec[:] = f(x)
    return convert(Cint, 1)
end

const func_c = cfunction(func_, Cint, (Ptr{Int32}, Ptr{Int32}, Ptr{Float64}, Ptr{Float64}, Ptr{Void}))

function jacobian_(m_::Ptr{Int32}, n_::Ptr{Int32}, x_::Ptr{Float64}, fjac_::Ptr{Float64}, jac_::Ptr{Void})
    m = unsafe_load(m_)
    n = unsafe_load(n_)
    x = pointer_to_array(x_, (n,))
    fjac = pointer_to_array(fjac_, (m,n))
    jac = unsafe_pointer_to_objref(jac_)::Function
    fjac[:,:] = jac(x)
    return convert(Cint, 1)
end

const jacobian_c = cfunction(jacobian_, Cint, (Ptr{Int32}, Ptr{Int32}, Ptr{Float64}, Ptr{Float64}, Ptr{Void}))
    
function avv_(m_::Ptr{Int32}, n_::Ptr{Int32}, x_::Ptr{Float64}, v_::Ptr{Float64}, acc_::Ptr{Float64}, avv_::Ptr{Void})
    m = unsafe_load(m_)
    n = unsafe_load(n_)
    x = pointer_to_array(x_, (n,))
    v = pointer_to_array(v_, (n,))
    acc = pointer_to_array(acc_, (m,))
    avv = unsafe_pointer_to_objref(avv__)::Function
    acc[:] = avv(x,v)
    return convert(Cint, 1)
end

const avv_c = cfunction(avv_, Cint, (Ptr{Int32}, Ptr{Int32}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Void}))

function callback_(m_::Ptr{Int32}, n_::Ptr{Int32}, x_::Ptr{Float64}, v_::Ptr{Float64}, a_::Ptr{Float64}, fvec_::Ptr{Float64}, fjac_::Ptr{Float64}, acc_::Ptr{Float64}, lam_::Ptr{Float64}, dtd_::Ptr{Float64}, fvec_new_::Ptr{Float64}, accepted_::Ptr{Int32}, info_::Ptr{Int32}, cb_::Ptr{Void})
    m = unsafe_load(m_)
    n = unsafe_load(n_)
    x = pointer_to_array(x_, (n,))
    v = pointer_to_array(v_, (n,))
    a = pointer_to_array(a_, (n,))
    fvec = pointer_to_array(fvec_, (m,))
    fjac = pointer_to_array(fvec_, (m,n))
    acc = pointer_to_array(acc_, (m,))
    lam = unsafe_load(lam_)
    dtd = pointer_to_array(dtd_, (n,n))
    fvec_new = pointer_to_array(fvec_new_, (m,))
    accepted = unsafe_load(accepted_)
    info = pointer_to_array(info_, (1,))
    cb = unsafe_pointer_to_objref(cb_)::Function
    info[1] = cb(m, n, x, v, a, fvec, fjac, acc, lam, dtd, fvec_new, accepted, info)
    return convert(Cint, 1)
end    

const callback_c = cfunction(callback_, Cint, (Ptr{Int32}, Ptr{Int32}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Int32}, Ptr{Int32}, Ptr{Void}))


function jacobian_dummy(x)
    return
end

function avv_dummy(x,v)
    return
end

function callback_dummy(m, n, x, v, a, fvec, fjac, acc, lam, dtd, fvec_new, accepted, info)
    return 0
end

## Wrapper around the main geodesiclm routine
function geodesiclm(func::Function, x::Array{Float64}, m::Integer, n::Integer;
                    jacobian::Function = jacobian_dummy, avv::Function = avv_dummy, callback::Function = callback_dummy,
                    center_diff::Bool = false, h1::Float64 = sqrt(eps()), h2::Float64 = eps()^0.25,
                    damp_mode::Integer = 1,
                    maxiter::Integer = 500, maxfev::Integer = 0, maxjev::Integer = 0, maxaev::Integer = 0, maxlam::Float64 = -1.0, minlam::Float64 = -1.0,
                    artol::Float64 = 1e-3, Cgoal::Float64 = sqrt(eps()), gtol::Float64 = sqrt(eps()), 
                    xtol::Float64 = sqrt(eps()), xrtol::Float64 = -1.0, ftol::Float64 = sqrt(eps()), frtol::Float64 = -1.0, 
                    print_level::Integer = 1, print_unit::Integer = 6,
                    imethod::Integer = 0, iaccel::Integer =1, ibold::Integer = 0, ibroyden::Integer = 0,
                    initialfactor::Float64 = 0.001, factoraccept::Float64 = 3.0, factorreject::Float64 = 2.0, avmax::Float64 = 0.75
                    )

    analytic_jac = (jacobian != jacobian_dummy)
    analytic_avv = (avv != avv_dummy)

    fvec = zeros(m)
    fjac = zeros(m,n)

    # we need to get the value of info back from geodesiclm.  This seems to be the easiest way to make the changed value accessible after the call
    info = [convert(Int32,0)]

    # at the moment, we don't pass a dtd array
    dtd = eye(n)
    
    niters = [convert(Int32, 0)]
    nfev = [convert(Int32,0)]
    njev = [convert(Int32,0)]
    naev = [convert(Int32, 0)]
    converged = [convert(Int32,0)]
    
    ccall( (:geodesiclm_wrapper, geodesiclmlib), Void, 
          (Ptr{Void}, Ptr{Void}, Ptr{Void}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Int32}, Ptr{Int32}, Ptr{Void}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Any, Any, Any, Any), 
    func_c, jacobian_c, avv_c, x, fvec, fjac, &convert(Int32, n), &convert(Int32, m), callback_c, info,
    &analytic_jac, &analytic_avv, &center_diff, &h1, &h2, dtd, &convert(Int32, damp_mode), 
    niters, nfev, njev, naev, 
    &convert(Int32, maxiter), &convert(Int32,maxfev), &convert(Int32,maxjev), &convert(Int32,maxaev),
    &maxlam, &minlam, &artol, &Cgoal, &gtol, &xtol, &xrtol, &ftol, &frtol, converged,
    &convert(Int32,print_level), &convert(Int32,print_unit),
    &convert(Int32,imethod), &convert(Int32,iaccel), &convert(Int32,ibold), &convert(Int32,ibroyden),
    &initialfactor, &factoraccept, &factorreject, &avmax, func, jacobian, avv, callback)

    return x, {"converged"=>converged[1], "iters"=>[niters[1], nfev[1], njev[1], naev[1]], "msg"=>converged_msg[converged[1]], "fvec"=>fvec, "fjac"=>fjac}
end

export geodesiclm

end # module

