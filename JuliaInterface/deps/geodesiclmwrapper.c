void geodesiclm_wrapper(
			void (*func)(int*,int*,float*,float*,void*), /* func with thunk argument */
			void (*jacobian)(int*, int*, float*, float*,void*), /* jacobian with thunk argument */
			void (*avv)(int*, int*, float*, float*, float*, void*), /* avv with thunk argument */
			float *x, float *fvec, float *fjac, int *n, int *m,
			void (*callback)(int*,int*,float*,float*,float*,float*,float*,float*,float*,float*,float*,int*,int*,void*),
			/* callback with thunk argument */
			int *info, int *analytic_jac, int *analytic_avv, int *center_diff, 
			float *h1, float *h2, float *dtd, int *damp_mode,
			int *niters, int *nfev, int *njev, int *naev,
			int *maxiter, int *maxfev, int *maxjev, int *maxaev, float *maxlam, float *minlam,
			float *artol, float *Cgoal, float *gtol, float *xtol, float *xrtol, float *ftol, float *frtol,
			int *converged, int *print_level, int *print_unit,
			int *imethod, int *iaccel, int *ibold, int *ibroyden,
			float *initialfactor, float *factoraccept, float *factorreject, float *avmax,
			/* thunk parameters */
			void *functhunk, void *jacobianthunk, void *avvthunk, void *callbackthunk
			){

  void newfunc(int *m, int *n, float *x, float *fvec){
    func(m, n, x, fvec, functhunk);
  };

  void newjacobian(int *m, int *n, float *x, float *fjac){
    jacobian(m, n, x, fjac, jacobianthunk);
  };
  
  void newavv(int *m, int *n, float *x, float *v, float *acc){
    avv(m, n, x, v, acc, avvthunk);
  };

  void newcallback(int *m, int *n, float*x, float*v, float *a, float *fvec, float *fjac, float *acc,
		   float *lam, float *dtd, float *fvec_new, int *accepted, int *info){
    callback(m,n,x,v,a,fvec,fjac,acc,lam,dtd,fvec_new,accepted,info,callbackthunk);
  };

  geodesiclm_(newfunc, newjacobian, newavv,
	      x, fvec, fjac, n, m, newcallback, info,
	      analytic_jac, analytic_avv, center_diff, h1, h2,
	      dtd, damp_mode,
	      niters, nfev, njev, naev, maxiter, maxfev, maxjev, maxaev, maxlam, minlam,
	      artol, Cgoal, gtol, xtol, xrtol, ftol, frtol,
	      converged, print_level, print_unit,
	      imethod, iaccel, ibold, ibroyden,
	      initialfactor, factoraccept, factorreject, avmax);
}
			
									
