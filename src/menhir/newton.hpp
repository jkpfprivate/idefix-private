// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef NEWTON_HPP_
#define NEWTON_HPP_

#include <petscsnes.h>
#include "menhir.hpp"
//#include "dataBlock.hpp"
//#include "problem.h"

#define SEPARATOR "---------------------------------------------------------------"

static PetscErrorCode SNESJacobianCallback(SNES, Vec, Mat, Mat, void*);
static PetscErrorCode SNESFunctionCallback(SNES, Vec, Vec, void *);
static PetscErrorCode SNESMonitorFunctionCallback(SNES, PetscInt, PetscReal, void*);

class Newton : public Menhir {
 public:
  Newton(int, char**);

  int Nconstraints;
  int Nmonitors;
  int Ntot;
  int Ntotloc;
  IdefixArray1D<real> constraints;
  IdefixArray1D<real> constraints0;
  IdefixArray1D<real> monitor;
  std::string  SNESRestartFile;
  std::string  SNESIterationFile;
  std::string  SNESResidualFile;
  std::string  SNESSolutionFile;
  std::string  action;
//  bool SNESMatrixFree;
  bool fixedPoint;

  PetscBool currentlyConstructingJacobian = PETSC_FALSE;

  /* Newton solver */
  PetscErrorCode Solve(void);
  
  /* Jacobian routine */
  PetscErrorCode SNESJacobian(SNES, Vec, Mat, Mat);
  
  /* Nonlinear function */
  PetscErrorCode SNESFunction(SNES, Vec, Vec);

  /* I/O and monitoring routines */
  PetscErrorCode SNESMonitorFunction(SNES, PetscInt, PetscReal);

  void print(IdefixArray4D<real> arr);
  IdefixArray4D<real> newtonGuess;

 private:
  
  /* Initial guess function */
  PetscErrorCode SNESInitialGuess(Vec, void*);
  
  void SNESReadVec(Vec, std::string);
  void SNESWriteVec(Vec, PetscInt, std::string);
  void ShowConfig(struct Scal*);
  
  void MapFieldForw(IdefixArray4D<real>, Vec);
  void MapFieldBack(Vec, IdefixArray4D<real>);
  void ProblemMonitor(IdefixArray4D<real>, IdefixArray1D<real>, IdefixArray1D<real>, int, int*);
  void PropagateField(real);
  void ComputeSNESFieldResidual(IdefixArray4D<real>, IdefixArray1D<real>, IdefixArray4D<real>);
  void ComputeSNESScalarResidual(IdefixArray4D<real>, IdefixArray1D<real>, IdefixArray1D<real>);

//  /* Constraints data to formulate the augmented system */
  IdefixArray4D<real> field;
  IdefixArray4D<real> field0;
  IdefixArray4D<real> field1;
  int newtonGuessIts;
  int nvar;
  
  /* Newton Flag */
  int NewtonFlag;
};
#endif //NEWTON_HPP_

