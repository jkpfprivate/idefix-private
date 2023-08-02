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

//using BodyForceFunc = void (*) (DataBlock &, const real t, IdefixArray4D<real>&);


class Forcing {
 public:
  Forcing(Input&, DataBlock*);
  void ComputeForcing(real);           ///< compute forcing field at current time t

  void ResetForcingTerm();            ///< fill the forcing field with zeros.

  void ShowConfig();                ///< Show the forcing configuration

  // Forcing term
  IdefixArray4D<real> forcingTerm;

//  // Whether we should skip gravity computation every n steps
//  int skipGravity{1};

  #if GEOMETRY == SPHERICAL
    real lmax;
    real mmax;
  #endif

 private:
//  friend class PlanetarySystem;

  DataBlock *data;

  real t_corr;
  real eps_Ylm;
  real eps_Slm;
  real eps_Tlm;
  OrnsteinUhlenbeckProcesses OUprocesses;

};

#endif // FORCING_FORCING_HPP_
