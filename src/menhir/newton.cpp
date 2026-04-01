// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#include <petscsnes.h>
//#include "util.h"
//#include "problem.h"
//#include "solvconf.h"
#include "idefix.hpp"
#include "dumpImage.hpp"
#include "dataBlock.hpp"
#include "newton.hpp"
#include "output.hpp"
#include "solvconf.h"

Newton::Newton(int argc, char* argv[]) : Menhir(argc, argv) {

  /* Variables required to compute function increments correctly */
  PetscBool currentlyConstructingJacobian = PETSC_FALSE;
  
  action = "Solve a nonlinear system using Newton method.";
  
  /* Default I/O */
  SNESRestartFile = "DNSfinal";
  SNESIterationFile = "Newtonit";
  SNESResidualFile = "Newtonres";
  SNESSolutionFile = "Newtonsol";
  
  /* Newton Flag */
  NewtonFlag = 0;

  fixedPoint = true;

  newtonGuessIts = 0; //????

  Nconstraints = 0;
  Nmonitors = 0;
  Ntotloc = data->np_tot[IDIR] * data->np_tot[JDIR] * data->np_tot[KDIR] * nvar + Nconstraints;
  Ntot = Ntotloc;

  /* Newton Guess - Can be used by problem routines to compute residuals
     for extra constraints variables (e. g. Viswanath's ortho condition) */
  fieldNewtonGuess = IdefixArray4D<real>("fieldNewtonGuess", nvar, data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);
  fieldResidual = IdefixArray4D<real>("fieldResidual", nvar, data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);
  fieldInitial = IdefixArray4D<real>("fieldInitial", nvar, data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);
  fieldMonitor = IdefixArray4D<real>("fieldMonitor", nvar, data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);
  fieldSaved = IdefixArray4D<real>("fieldSaved", nvar, data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);
  constraints = IdefixArray1D<real>("constraints", Nconstraints);
  constraints0 = IdefixArray1D<real>("constraints0", Nconstraints);
  monitor = IdefixArray1D<real>("monitor", Nmonitors);

}

/* ------------------------------------------------------------------------ 
                         Nonlinear Newton Solver
   ------------------------------------------------------------------------ */

PetscErrorCode Newton::Solve(){

  /* Variables: 
     snes        - nonlinear solver
     ksp         - Krylov method (ksptolerance = method tolerance)
     pc          - preconditioner
     X,F         - X is the state vector, F serves to calculate function
     J           - Jacobian matrix
     Note: even the matrix-free implementation requires
     declaring and initializing a J matrix formally
     its         - # iterations for convergence
     ierr        - error code */
  
  SNES             snes;
  SNESLineSearch   linesearch;
  KSP              ksp;
  PC               pc;
  Vec              X,F;
  Mat              J;
  PetscInt         i,its;
  PetscErrorCode   ierr;
  SNESConvergedReason convReason;
  char solutionf[53];

  /* Set Newton flag to 1 */
  NewtonFlag = 1;

  /* There is no implementation of the explicit jacobian version with mpi ! */
  if ((idfx::psize > 1) && (!SNESMatrixFree)){
    SETERRQ(MPI_COMM_WORLD,0,"No MPI Version of the Newton solver with explicit jacobian construction !\n");
  }

  /* Create nonlinear solver context */
  ierr = SNESCreate(MPI_COMM_WORLD,&snes);

  /* Nonlinear solver configuration */
  ierr = SNESSetType(snes,SNESSolverType);
  ierr = SNESGetLineSearch(snes,&linesearch); 
  ierr = SNESLineSearchSetType(linesearch,SNESLinesearch);
//  ierr = SNESLineSearchSet(snes,SNESLinesearch,PETSC_NULLPTR); 
//  ierr = SNESLineSearchSetParams(snes,SNESLSalpha,SNESLSmaxstep);
  ierr = SNESSetTolerances(snes,SNESatol,SNESrtol,SNESstol,SNESNitmax,SNESfevalmax); 
  ierr = SNESMonitorSet(snes,SNESMonitorFunctionCallback,this,PETSC_NULLPTR); 

  PetscViewerAndFormat *vf;
  PetscViewerAndFormatCreate(PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_DEFAULT, &vf);

  /* Krylov Method configuration */
  ierr = SNESGetKSP(snes,&ksp); 
  ierr = KSPSetType(ksp,SNESKSPMethod);  
  ierr = KSPSetTolerances(ksp,SNESKSPrtol,SNESKSPatol,SNESKSPdtol,SNESKSPNitmax);  
  ierr = KSPMonitorSet(ksp,(PetscErrorCode (*)(KSP, PetscInt, PetscReal, void*)) SNESKSPMonitor,vf,PETSC_NULLPTR);
  ierr = KSPGMRESSetRestart(ksp,SNESKSPNrestart);

  /* KSP Preconditioner configuration */
  ierr = KSPGetPC(ksp,&pc);  
  ierr = PCSetType(pc,SNESKSPPCMethod);  

  /* Create vector data structures; set corresponding routines */
  ierr = VecCreateMPI(MPI_COMM_WORLD,Ntotloc,Ntot,&X); 
  ierr = VecDuplicate(X,&F) ; 

  /* Set function evaluation routine and vector */
  ierr = SNESSetFunction(snes,F,SNESFunctionCallback,this);

  /* Create Jacobian Matrix  */ 
  if (SNESMatrixFree){
    ierr = MatCreateSNESMF(snes,&J); 
    ierr = MatMFFDSetType(J,SNESMFIncrement);  
    ierr = MatMFFDSetFunctionError(J,SNESrerror);  
    ierr = MatMFFDDSSetUmin(J,SNESXmin);  
  }
  else{ierr = MatCreateDense(MPI_COMM_WORLD,Ntot,Ntot,Ntot,Ntot,PETSC_NULLPTR,&J);}

  /* Set Jacobian evaluation routine */
  /* SNESFunction is used to estimate [J]X in the matrix-free version*/
  ierr = SNESSetJacobian(snes,J,J,SNESJacobianCallback,this);

  /* Generate initial guess */
  ierr = SNESInitialGuess(X,PETSC_NULLPTR);

  /* Solve problem */
  ierr = SNESSolve(snes,PETSC_NULLPTR,X); 

  /* Final I/O */  
  ierr = SNESGetConvergedReason(snes,&convReason);
  if (convReason > 0){
    printf("Newton Solver has converged !\n");
//    MPI_Printf("Newton Solver has converged !\n");
//    sprintf(solutionf,"%s",SNESSolutionFile);
    SNESWriteVec(X,newtonGuessIts,SNESSolutionFile);
  }
  else{printf("Warning: Newton Solver has not converged !\n");}
              
  /* ierr = SNESGetIterationNumber(snes,&its);  */
  printf("Total number of Newton iterations = %3d\n",
	     newtonGuessIts-1);

  /*  Free work space */
  ierr = MatDestroy(&J); ierr = VecDestroy(&X); 
  ierr = VecDestroy(&F); ierr = SNESDestroy(&snes); 

  return(0);
}

/* ------------------------------------------------------------------------ 
                              Form Initial Guess
   ------------------------------------------------------------------------ */

PetscErrorCode Newton::SNESInitialGuess(Vec X,void *ctx){

  /* Read State Vector from a file */
  data->dump->Read(*output, 0);
//  ReadField(field,restartf);
//  ReadScal(&scal,restartf);

//  mysetup->InitFlow(*data);
  data->SetBoundaries();
  data->Validate();
  Copy(fieldInitial,data->hydro->Vc);

  /* its=-1 means that Newton has been called by another action i.e. continuation.
     We don't want to perform the next two operations in that case - They should 
     have been performed earlier */

  /* If we are not dealing with a restart or a continuation, obtain guesses for 
     the constraints variables from problem.c */
//  if (scal.its == 0){GetCustomScalars(&scal);}

  /* Change problem parameters according to the information just read */
//  SetProblemParameters(&scal);

//  if (scal.its >= 0){
//    /* Compute the scaling Factors for the mappings - This is mandatory ! */
//    ComputeRescalingFactors(&field,&scal);
//  }
//  if (scal.its == -1){scal.its = 0;}
  
  /* Once the scaling factors have been calculated, we can map fields and 
     scalars into X (these factors are used in the mapping ! */
#if MHD == YES
  data->hydro->boundary->ReconstructVcField(data->hydro->Vc);
#endif //MHD == YES
//  ShowVc();
  MapFieldForw(data->hydro->Vc,X); //MapScalForw(&scal,X);
//  ShowX(X);

  /* Initialize the constraints with this initial guess */
  MapFieldBack(X,fieldNewtonGuess);
//  MapScalBack(X,&newtonguess.scal);
//  newtonguess.scal.its=scal.its;

  /* Write Newton action information to file */
//  ShowConfig(&scal);
  
  return(0);
}

/* ------------------------------------------------------------------------ 
                          Compute nonlinear function
   ------------------------------------------------------------------------ */
PetscErrorCode Newton::SNESFunction(SNES snes,Vec X,Vec F) {
//idfx::cout << "Entering SNESFunction" << std::endl;
  /* Variables */  
  KSP ksp;
  PetscErrorCode ierr;
  double *crhs,*extravars;
  PetscInt i,its;
  KSPConvergedReason convReason;

  /* Book-keeping of running Newton guess after the Krylov solver calls
     to the function: this is required to compute the correct residual before
     the next newton iteration */
  ierr = SNESGetKSP(snes,&ksp); 
  ierr = KSPGetConvergedReason(ksp,&convReason); 

//  if (true){    
  if ((convReason > 0) && (!currentlyConstructingJacobian)){    
    MapFieldBack(X,fieldNewtonGuess);
//    MapScalBack(X,&newtonguess.scal);
  }

  /* Calculate function - start by transforming back to field/scalar structures */
  MapFieldBack(X,fieldInitial); //MapScalBack(X,&scal0);
  
  /* Compute Field components of the Residual of the total system 
     i.e. integrate in time, take difference etc. */
  ComputeSNESFieldResidual();

//  /* ! Solvability conditions (i.e. Viswanath-like condition etc.) */
//  ComputeSNESScalarResidual(fieldInitial,constraints0,constraints);

  /* Get back to a state vector representation */
#if MHD == YES
  data->hydro->boundary->ReconstructVcField(data->hydro->Vc);
#endif //MHD == YES
  MapFieldForw(fieldResidual,F); //MapScalForw(&scal,F);

//idfx::cout << "Exiting SNESFunction" << std::endl;
  return 0;
}

static PetscErrorCode SNESFunctionCallback(SNES snes, Vec x, Vec F, void *ctx)
{
    auto *self = static_cast<Newton*>(ctx);
    return self->SNESFunction(snes, x, F);
}

///* ------------------------------------------------------------------------ 
//                         Evaluate Jacobian matrix
//   ------------------------------------------------------------------------ */
//
PetscErrorCode Newton::SNESJacobian(SNES snes,Vec X,Mat J, Mat B) {
//idfx::cout << "Entering SNESJacobian" << std::endl;
  /* Input/Output variables */
  PetscErrorCode ierr;
  PetscScalar *X_v;

  /* extra variables for FD estimate of Jacobian */
  int i,k;
  PetscInt *idx;
  PetscScalar *A;
  Vec Xinc,F,Finc;
  PetscScalar *Xinc_v,*F_v,*Finc_v;

  /* FD estimate of Jacobian - coded here just for comparison with Carlo */

  if (!SNESMatrixFree){

    idx = (PetscInt *) malloc( sizeof(PetscInt) * Ntot);
    A = (PetscScalar *) malloc( sizeof(PetscScalar) * Ntot * Ntot) ;
    for (i=0 ; i < Ntot ; i++){
      idx[i] = 0. ;
      for (k=0 ; k < Ntot ; k++){
	A[k+Ntot*i]=0. ;
      }
    }
      
    currentlyConstructingJacobian=PETSC_TRUE;
    ierr = VecDuplicate(X,&Xinc); 
    ierr = VecDuplicate(X,&F); 
    ierr = VecDuplicate(X,&Finc); 
    
    ierr = SNESFunction(snes,X,F);
    for (k=0 ; k < Ntot ; k++){
      ierr = VecGetArray(X,&X_v); 
      ierr = VecGetArray(Xinc,&Xinc_v); 
      for (i=0 ; i < Ntot ; i++){ Xinc_v[i]=X_v[i];}

     /* small increment */
      if (X_v[k] != 0.){Xinc_v[k]=1.001*X_v[k];}
      else {Xinc_v[k]=0.0001;}

      ierr = VecRestoreArray(Xinc,&Xinc_v); 
      ierr = SNESFunction(snes,Xinc,Finc); 
      ierr = VecGetArray(Xinc,&Xinc_v); 
      ierr = VecGetArray(F,&F_v); 
      ierr = VecGetArray(Finc,&Finc_v); 
      
      /* WARNING: do we really want row-major format for Aik ? */
      for (i=0 ;  i < Ntot ; i++){A[k+Ntot*i]=(Finc_v[i]-F_v[i])/(Xinc_v[k]-X_v[k]);}

      ierr = VecRestoreArray(X,&X_v); 
      ierr = VecRestoreArray(Xinc,&Xinc_v); 
      ierr = VecRestoreArray(F,&F_v); 
      ierr = VecRestoreArray(Finc,&Finc_v); 
    }

    for (i=0 ; i < Ntot ; i++){idx[i]=i;}
    ierr = MatSetValues(J,Ntot,idx,Ntot,idx,A,INSERT_VALUES); 
    currentlyConstructingJacobian=PETSC_FALSE;
    free (A) ; free(idx) ;

  }

  /* END OF FD JACOBIAN SECTION */

  /* Assemble matrix - Warning: these calls are required independently of
     the choice of a Matrix-Free implementation !!! */

  ierr = MatAssemblyBegin(J,MAT_FINAL_ASSEMBLY); 
  ierr = MatAssemblyEnd(J,MAT_FINAL_ASSEMBLY); 

//idfx::cout << "Exiting SNESJacobian" << std::endl;
  return 0;
}

static PetscErrorCode SNESJacobianCallback(SNES snes, Vec x, Mat J, Mat P, void *ctx)
{
    auto *self = static_cast<Newton*>(ctx);
    return self->SNESJacobian(snes, x, J, P);
}

/* ------------------------------------------------------------------------ 
                         Monitor the Newton Solver
   ------------------------------------------------------------------------ */

static PetscErrorCode SNESMonitorFunctionCallback(SNES snes, PetscInt its, PetscReal rnorm, void *ctx)
{
    auto *self = static_cast<Newton*>(ctx);
    return self->SNESMonitorFunction(snes, its, rnorm);
}

PetscErrorCode Newton::SNESMonitorFunction(SNES snes, PetscInt its, PetscReal rnorm){

  /* on input, rnorm this is the residual norm apparently */
  PetscErrorCode ierr;
  PetscScalar *F_v;
  Vec X,F;
  double xnorm;
  PetscInt i,it;
  int Nmonitors;
//  struct Field field;
//  struct Scal scal,monitor;
  KSP ksp ;
  char iterationf[50], residualf[50];
//  std::string iterationf, residualf;
  FILE *ht;

//  AllocateField(&field); AllocateScal(&scal); AllocateScal(&monitor);

  /* Use the actual total number of iteration - Different from its if
     we have made a restart from a previous unfinished snes solve */
  it = newtonGuessIts;

  /* get current Newton guess and compute its norm */
  ierr = SNESGetSolution(snes,&X); 

  VecNorm(X, NORM_2, &xnorm);
 
  /* get current residual */
  ierr = SNESGetFunction(snes,&F,PETSC_NULLPTR,PETSC_NULLPTR);

  /* User decides what he wants to output to terminal */
  MapFieldBack(X,fieldMonitor); //MapScalBack(X,&scal);
  ProblemMonitor(fieldMonitor, constraints, monitor, Nconstraints, &Nmonitors);
  
  /* Check the residual norm and state vector */
  idfx::cout << std::endl;
  idfx::cout << "It =" << std::setw(4) << it << " || Residual Norm = " << rnorm << std::endl;
  idfx::cout << "         ||    Guess Norm = " << xnorm << std::endl;
  idfx::cout << "         ||  State Vector = ";
  for (i=0 ; i < Nmonitors ; i++){
    idfx::cout << monitor(i);
  }
  idfx::cout << std::endl;

  if (idfx::prank==0){
    if (it == 0){ht=fopen("Newtonconvergence.dat","w");}
    else{ht=fopen("Newtonconvergence.dat","a");}
    fprintf(ht,"%3d %20.16f",it,rnorm);
//    MPI_Fprintf(ht,"%3d %20.16f",it,rnorm);
    for (i=0 ; i < Nmonitors ; i ++){
//      MPI_Fprintf(ht," %20.16f",monitor(i));
      fprintf(ht," %20.16f",monitor(i));
    }
//    MPI_Fprintf(ht,"\n");
    fprintf(ht,"\n");
    fclose(ht);
  }

  /* Dump Current Newton guess */
//  sprintf(iterationf,"%s%3.3d",SNESIterationFile,it);
//  SNESWriteVec(X,it,iterationf);
  SNESWriteVec(X,it,SNESIterationFile);

  /* Dump Current Residual to file */
//  sprintf(residualf,"%s%3.3d",SNESResidualFile,it);
//  SNESWriteVec(F,it,residualf);
  SNESWriteVec(F,it,SNESResidualFile);

  /* Increment the current Newton guess iteration number by one */
  newtonGuessIts += 1;

  return 0;
}

///* ------------------------------------------------------------------------ 
//                       Write State Vectors to file
//              They contain physical fields + extra variables
//   ------------------------------------------------------------------------ */
//void Newton::SNESReadVec(Vec X, std::string myfile){
//
//  DumpImage image(myfile, &(*data));
//  MapFieldForw(data->hydro->Vc, X);
////  MapScalForw(data, X);
//  //destroy data?
//}

void Newton::SNESWriteVec(Vec X, PetscInt it, std::string myFile){
  Copy(fieldSaved,data->hydro->Vc);
//ShowVc();
  MapFieldBack(X, data->hydro->Vc);
#if MHD == YES
  data->hydro->boundary->ReconstructVcField(data->hydro->Vc);
#endif //MHD == YES
  data->dump->Write(*output, myFile, it);
  Copy(data->hydro->Vc, fieldSaved);
//  MapScalBack(X,&scal); 
//  scal.its=it; /* set scalars, including current snes iteration */
//  scal.time=0.;
//  WriteScal(&scal,myfile);
}

/* ------------------------------------------------------------------------ 
                  Dump Action Information to terminal & file
   ------------------------------------------------------------------------ */
void Newton::ShowConfig(struct Scal *scal){

  FILE *ht;
  int i;

  idfx::cout << SEPARATOR << std::endl;
  idfx::cout << action << std::endl;
  idfx::cout << "Initial guess in file " << SNESRestartFile << std::endl;
  idfx::cout << "Extra variables guesses:" << std::endl;
  if (Nconstraints !=0){
//    for (i=0;i<Nconstraints; i++){idfx::cout << scal->vars[i];}
    idfx::cout << std::endl;
  }
  idfx::cout << "Solution saved in file" << SNESSolutionFile << std::endl;
  idfx::cout << SEPARATOR << std::endl;
}

//void Newton::ShowX(Vec X) {
//  int nz, ny, nx;
//  nz = data->np_tot[KDIR];
//  ny = data->np_tot[JDIR];
//  nx = data->np_tot[IDIR];
//  const PetscScalar *X_v;
//  int ierr;
//  ierr = VecGetArrayRead(X,&X_v);
//
//  idefix_for("MapFieldBack", 0,nvar,
//             0,data->np_tot[KDIR],
//             0,data->np_tot[JDIR],
//             0,data->np_tot[IDIR],
//    KOKKOS_LAMBDA (int v, int k, int j, int i) {
//      idfx::cout << X_v[UNRAVEL(v,k,j,i)] << std::endl;
//    });
//  ierr = VecRestoreArrayRead(X,&X_v);
//}
