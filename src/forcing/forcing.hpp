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

#if GEOMETRY == SPHERICAL
  #include "vsh.hpp"
#endif
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


 private:
//  friend class PlanetarySystem;

  DataBlock *data;

  real t_corr;
  real frms_Ylm;
  real frms_Slm;
  real frms_Tlm;
  OrnsteinUhlenbeckProcesses OUprocesses;

  #if GEOMETRY == SPHERICAL
//    Vsh vsh;
    std::unique_ptr<Vsh> vsh;
    real lmax;
    real mmax;
  #endif

};

#endif // FORCING_FORCING_HPP_
