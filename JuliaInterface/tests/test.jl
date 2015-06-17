import GeodesicLM


# Rosenbrock test model
function r(x)
    return [1.0 - x[1], 10.0*(x[2] - x[1]^2)]
end

function j(x)
    # return [ -1.0 -20.0*x[1]; 0.0 10.0]
    return [ -1.0 0.0; -20.0*x[1] 10.0]
end

function avv(x, v)
    return [0.0, -20.0*v[1]^2]
end

## Call the wrapper and get informtion about the fit

x, info = GeodesicLM.geodesiclm(r, [-2.0, 4.0], 2, 2)
println("x = $x")
println("""msg = $(info["msg"])""")
fvec = info["fvec"]
println("Final Cost = $(sum(fvec.*fvec)/2)")

#= Result of running this script:

 Optimizing with Geodesic-Levenberg-Marquardt algorithm, version 1.0
 Method Details:
   Update method:              0
   acceleration:               1
   Bold method:                0
   Broyden updates:            0
   Initial Cost:       4.5000000000000000
 Optimization finished
 Results:
   Converged:    Cgoal reached              2
   Final Cost:    1.0899309828919105E-009
   Cost/DOF:                   Infinity
   niters:               13
   nfev:                 14
   njev:                 13
   naev:                 13
x = [0.9999533133323716,0.9999065822577443]
Final Cost = 1.0899309828919105e-9

=#
