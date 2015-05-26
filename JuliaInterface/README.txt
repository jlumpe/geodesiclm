README.txt for the Julia wrapper to the geodesic Levenberg-Marquardt routine.

To begin, compile the geodesicLM FORTRAN routine AND LINK IT to the lapack and blas libraries.  You should use the -shared and -fPIC options when compiling.  The current make file does not link, so you will need to do this by hand.  You can this, for example, with this command in the geodesicLM folder:

gfortran -shared -fPIC *.f *.f90 -o libgeodesiclm.so -LPATH_TO_LIBRARIES -llapack -lblas

where you replace PATH_TO_LIBRARIES with the appropriate path.

Move the resulting libgeodesiclm.so to the same folder with the julia wrapper.

GeodesicLM.jl provides one wrapper function to the FORTRAN routine.  It can be called minimally as

geodesiclm(func, x, m, n)

where func is a function with the signature:
func(m, n, x, fvec)
m = integer, size of fvec
n = integer, size of x
x = array of Float64 (size n), point at which the function is evaluated
fvec = array of Float 64 (size m), that should contain the values of the function at point x upon return
func should return a Cint

x is an array of Float64 values that are the initial guess for the fit
m = integer, number of residuals in the sum of squares
n = integer, number of parameters to fit

Other optional arguments can be supplied.  Most importantly:
jacobian: function to calculate the Jacobian matrix with the following signature
jacobian(m, n, x, fjac)
parameters are the same as func except fjac contains the jacobian upon exit
jacobian should return a Cint

Avv: Function to calculate the directional second derivative with the following signature
Avv(m,n,x,v,acc)
m,n,x are same as func
v = direction of directional derivative
acc contains the second directional derivative upon exit
Avv should return a Cint

Most other keyword arguments are the same as in the FORTRAN documentation

Upon exit, the routine returns the final parameter estimates and a dictionary with information about the fit.

