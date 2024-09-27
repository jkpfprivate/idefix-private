// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef FORCING_FORCING_HPP_
#define FORCING_FORCING_HPP_

#include "idefix.hpp"
#include "input.hpp"
#include "OrnsteinUhlenbeckProcess.hpp"

class DataBlock;

class Forcing {
 public:
  Forcing(Input&, DataBlock*);
  void ComputeForcing(real);           ///< compute forcing field at current time t

  void ResetForcingTerm();            ///< fill the forcing field with zeros.

  void ShowConfig();                ///< Show the forcing configuration

  // Forcing modes
  IdefixArray4D<real> forcingModes;
  // Forcing term
  IdefixArray4D<real> forcingTerm;
  OrnsteinUhlenbeckProcesses OUprocesses;

//  // Whether we should skip gravity computation every n steps
//  int skipGravity{1};

  int write;
 private:

  DataBlock *data;
  int seed;

  real t_corr;
  #if GEOMETRY == SPHERICAL & VSH == YES
    real eps_Ylm;
    real eps_Slm;
    real eps_Tlm;
  #endif

};

#endif // FORCING_FORCING_HPP_
