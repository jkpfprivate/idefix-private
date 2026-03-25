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

class Newton : public Menhir {
 public:
  Newton(int, char**);

  int Nconstraints;
  int Nmonitors;
  int Ntot;
  int Ntotloc;
  IdefixArray1D<real> constraints;
  IdefixArray1D<real> monitor;
  std::string  SNESRestartFile;
  std::string  SNESIterationFile;
  std::string  SNESResidualFile;
  std::string  SNESSolutionFile;
  std::string  action;
//  bool SNESMatrixFree;
  bool fixedPoint;

  PetscBool CurrentlyConstructingJacobian = PETSC_FALSE;
  private:
  /* Newton solver */
  PetscErrorCode Solve(void);
  
  /* Initial guess function */
  PetscErrorCode SNESInitialGuess(Vec, void*);
  
  /* Nonlinear function */
  static PetscErrorCode SNESFunction(SNES, Vec, Vec, void*);
  
  /* Jacobian routine */
  static PetscErrorCode SNESJacobian(SNES, Vec, Mat, Mat, void*);
  
  /* I/O and monitoring routines */
  static PetscErrorCode SNESMonitorFunc(SNES, PetscInt, PetscReal, void*);
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
  IdefixArray4D<real> newtonGuess;
  int newtonGuessIts;
  int nvar;
  
  /* Newton Flag */
  int NewtonFlag;
};
#endif //NEWTON_HPP_

