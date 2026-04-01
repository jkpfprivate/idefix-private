/* ------------------------ CONTINUATION SET-UP --------------------------- */

#define CONTNstep 5        /* Number of points requested on the continuation curve */
#define CONTvar .005 /* relative variation (.01=1%) of the user-defined 
			   continuation parameter in the first continuation step */
#define CONTrelax 1.       /* Relaxation factor for the arclength increment */

#define CONTrestart 0  /* 1 = restart from StartFile and corresponding TangentFile 
			  using the same continuation parameter as before */

/* ------------------------ END OF CONTINUATION SET-UP -------------------- */

/* --------------------------- NEWTON SOLVER SET-UP ----------------------- */

/* SNES Solver Configuration */

#define SNESSolverType SNESNEWTONLS /* or  SNESTR */
#define SNESLinesearch SNESLINESEARCHNONE /* or SNESLineSearchNo */
//#define SNESLinesearch SNESLineSearchCubic /* or SNESLineSearchNo */ 
#define SNESLSalpha .01 /* Parameters used by SNESLineSearchCubic...*/
#define SNESLSmaxstep 1e8

#define SNESrtol (double) 1e-8
#define SNESatol (double) 1e-10 /* or PETSC_DEFAULT */
#define SNESstol (double) 1e-20 /* or PETSC_DEFAULT */
#define SNESNitmax 100
#define SNESfevalmax 100000 

#define SNESMonitorFunction SNESMonitorFunction  /* or SnesMonitorDefault */

/* KSP Solver Configuration */

#define SNESKSPMethod KSPGMRES /* for methods see petsksp.h */
#define SNESKSPMonitor KSPMonitorResidual
#define SNESKSPPCMethod PCNONE

/* rtol is always what sets the convergence of KSP, for atol <<<<1 
/* Important note: rtol is a criterion on norm(r_i)/norm(r_0), 
/* where r_i is the residual after the i-th Krylov iteration 
/* and norm(r_0) is simply the Newton residual (rhs) in the  
/* linear Jacobian system  */

#define SNESKSPrtol (double) 1e-5
#define SNESKSPatol (double) 1e-50 
#define SNESKSPdtol (double) 100.0 
#define SNESKSPNrestart Ntot
#define SNESKSPNitmax 50

/* Matrix-Free Solver Configuration */

#define SNESMatrixFree 1 /* Set it to 0 to build Jacobian explicitly */

/* MATMFFD_DS: Dennis and Schnabel method: same as channelflow.org  
/* or MATMFFD_WP for Walker and Pernice Method (See Petsc manual p.98) 
/* Walker - Pernice seems to give a much more "stable" increment */

#define SNESMFIncrement MATMFFD_DS /* or MATMMFD_WP */

/* IMPORTANT WARNING: the parameter rerror sets the magnitude of  
/* the increment vector in the calculation of action of the Jacobian  
/* on vectors (see Petsc manual p. 98). The increment is ~rerror*a,
/* where a is the rhs of the linear jacobian system to solve.
/* Having a large enough rerror here is required to avoid issues with 
/* machine precision when estimating the jacobian action by some form 
/* of finite difference approximation. Having a small enough value of 
/* rerror is required to get a good FD approximation, though.
/* Empirically 1.e-4 to 1.e-3 works well (petsc default = 1e-7) */

#define SNESrerror (double) 1e-5

/* Xmin - minimum order of mag. of X_i under which Xmin is used instead 
   of X for the calculation of the jacobian increment h  */

#define SNESXmin (double) 1e-9

/* --------------------- END OF NEWTON SOLVER SET-UP ---------------------- */

/* ------------------- STABILITY SOLVER CONFIGURATION --------------------- */

/* Perturbations */

#define EPSpertmode  0 /* Choose eigenmode to perturb Newton state with */
#define EPSpertamp  -0.001 /* Perturbation relative amplitude */

/* Eigenvalue Solver Configuration */

#define EPSProblem EPS_NHEP  /* Non-Hermitian or Hermitian Matrix */
#define EPSMethod EPSARNOLDI /* or EPSKRYLOVSCHUR */

#define EPStol 1e-10
#define EPSitmax 20

#define EPSNeigen 8 /* Number of requested eigenmodes (the solver will probably give you more of them) */

#define EPSUseST  0 /* Use spectral transformation */
#define EPSSTType STSHIFT
#define EPSSHIFT -6.5946187419741804 
#define EPSWhich EPS_LARGEST_MAGNITUDE  /* how eigenvalues should be looked for */

/* Matrix-Free Solver Configuration */

#define EPSMatrixFree 1 /* Set it to 0 to build Jacobian explicitly */

/* MATMFFD_DS: Dennis and Schnabel method: same as channelflow.org  
/* or MATMFFD_WP for Walker and Pernice Method (See Petsc manual p.98) 
/* Walker - Pernice seems to give a much more "stable" increment */

#define EPSMFIncrement MATMFFD_WP /* or MATMMFD_WP */

/* IMPORTANT WARNING: the parameter rerror sets the magnitude of  
/* the increment vector in the calculation of action of the Jacobian  
/* on vectors (see Petsc manual p. 98). The increment is ~rerror*a,
/* where a is the rhs of the linear jacobian system to solve.
/* Having a large enough rerror here is required to avoid issues with 
/* machine precision when estimating the jacobian action by some form 
/* of finite difference approximation. Having a small enough value of 
/* rerror is required to get a good FD approximation, though. */

#define EPSrerror (double) 1e-8

/* Xmin - minimum order of mag. of X_i under which Xmin is used instead 
   of X for the calculation of the jacobian increment h  */

#define EPSXmin   (double) 1e-20

/* ---------------- END OF STABILITY SOLVER CONFIGURATION ----------------- */
