#include <petscsnes.h>
//#include "util.h"
//#include "problem.h"
//#include "solvconf.h"
#include "idefix.hpp"
#include "dumpImage.hpp"
#include "dataBlock.hpp"
#include "newton.hpp"
#include "output.hpp"

#define UNRAVEL(v,k,j,i) v*nvar*nz*ny*nx + k*nz*ny*nx + j*ny*nx + i*nx

void Newton::MapFieldForw(IdefixArray4D<real> Vc, Vec X) {
//  int nz, ny, nx;
//  nz = data->np_int[KDIR];
//  ny = data->np_int[JDIR];
//  nx = data->np_int[IDIR];
//  const PetscScalar *X_v;
//  PetscErrorCode ierr;
//  ierr = VecGetArrayRead(X,&X_v);
//
//  idefix_for("ConsToPrim", 0,nvar,
//             0,data->np_tot[KDIR],
//             0,data->np_tot[JDIR],
//             0,data->np_tot[IDIR],
//    KOKKOS_LAMBDA (int v, int k, int j, int i) {
//      Vc(v,k,j,i) = X_v[UNRAVEL(v,k,j,i)];
//    });
//
//  ierr = VecRestoreArrayRead(X,&X_v);
}

void Newton::MapFieldBack(Vec X, IdefixArray4D<real> Vc) {
//  int nz, ny, nx;
//  nz = data->np_int[KDIR];
//  ny = data->np_int[JDIR];
//  nx = data->np_int[IDIR];
//  PetscScalar *X_v;
//  PetscErrorCode ierr;
//  ierr = VecGetArray(X,&X_v);
//
//  idefix_for("ConsToPrim", 0,nvar,
//             0,data->np_tot[KDIR],
//             0,data->np_tot[JDIR],
//             0,data->np_tot[IDIR],
//    KOKKOS_LAMBDA (int v, int k, int j, int i) {
//      X_v[UNRAVEL(v,k,j,i)] = Vc(v,k,j,i);
//    });
//
//  ierr = VecRestoreArray(X,&X_v);
}

/* ------------------------------------------------------------------------ 
               Mandatory monitoring routines called at a higher level
   ------------------------------------------------------------------------ */

/* Provide a few monitoring quantities to higher level solvers  (Newton solver) */ 
void Newton::ProblemMonitor(IdefixArray4D<real> field, IdefixArray1D<real> constraints,
                    IdefixArray1D<real> monitor, int Nconstraints, int *Nmonitors){
  int i;
  /* We want to monitor Nmonitor scalar quantities */
//  *Nmonitors=Nconstraints+2;
  *Nmonitors=Nconstraints;

  /* Display cycle period, phase speed, continuation parameter */
  for (i=0 ; i<Nconstraints ; i++){monitor(i)=constraints(i);}

//  /* Display the first 2 field components */
//  monitor(Nconstraints)=field->vx[1];
//  monitor(Nconstraints+1)=field->bx[1];
}

void Newton::PropagateField(real tmax) {
  data->t = 0.;
  PerformDNS(tmax);
}

/* Calculate the field components of the rhs residual (Newton solver) */
void Newton::ComputeSNESFieldResidual(IdefixArray4D<real> field0, IdefixArray1D<real> constraints0, IdefixArray4D<real> field){

  int i;
  real dt,T,TLC,Cz,sz;
  IdefixArray4D<real> Tfield;
  
  if (fixedPoint){
    T=1.; /* integrate for a fixed time */
    /* Translate the result according to the phase speeds of the solution */
    Cz = constraints0(0);
  }
  else{
    /* Propagate fields to time T */
    T = constraints0[0];
    /* Translate the result according to the phase speeds of the solution */
    Cz = constraints0(1);
  }

//  SetProblemParameters(scal0); /* Update problem parameters and reset snoopy according to the input scalars */
  PropagateField(T);

//  sz = -Cz*T ; Translate(field,&Tfield,0.,0.,sz) ;

  /* Finally, compute the residual by taking the difference between
     V(T) translated back in space and V(0) */

  idefix_for("ConsToPrim",
             0,nvar,
             0,data->np_tot[KDIR],
             0,data->np_tot[JDIR],
             0,data->np_tot[IDIR],
    KOKKOS_LAMBDA (int v, int k, int j, int i) {
      field(v,k,j,i) = data->hydro->Vc(v,k,j,i) - field0(v,k,j,i);
    });
}


void Newton::ComputeSNESScalarResidual(IdefixArray4D<real> field0, IdefixArray1D<real> constraints0, IdefixArray1D<real> constraints){
}

