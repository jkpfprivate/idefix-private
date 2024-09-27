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

enum ForcingType {iso3D, iso2D, vsh, userDef};

class DataBlock;

class Forcing {
 public:
  Forcing(Input&, DataBlock*);

  void InitForcingModes();          ///< init forcing modes given its type

  void ComputeForcing(real);        ///< compute forcing field at current time t

  void ResetForcingTerm();          ///< fill the forcing field with zeros.

  void ShowConfig();                ///< Show the forcing configuration

  // Forcing modes
//  IdefixArray4D<Kokkos::complex<real>> forcingModes;
  IdefixArray4D<Kokkos::complex<real>> forcingModes;
  IdefixArray4D<real> forcingModesIdir;
  #if COMPONENTS >= 2
    IdefixArray4D<real> forcingModesJdir;
    #if COMPONENTS == 3
      IdefixArray4D<real> forcingModesKdir;
    #endif //COMPONENTS == 3
  #endif //COMPONENTS >= 2
  // Forcing term
  IdefixArray4D<real> forcingTerm;
  OrnsteinUhlenbeckProcesses OUprocesses;

//  // Whether we should skip gravity computation every n steps
//  int skipGravity{1};

  int write;
 private:

  DataBlock *data;
  int seed;

  int nForcingModes;
  int nForcingModesIdir;
  int nForcingModesJdir;
  int nForcingModesKdir;
  IdefixArray1D<real> tcorrs;
  IdefixArray1D<real> means;
  IdefixArray1D<real> epsilons;
//  #if GEOMETRY == SPHERICAL & VSH == YES
//    real eps_Ylm;
//    real eps_Slm;
//    real eps_Tlm;
//  #endif

////  #if GEOMETRY == CARTESIAN
//  bool have3Diso;
////  #endif //GEOMETRY == CARTESIAN
////  #if GEOMETRY == CARTESIAN or GEOMETRY == SPHERICAL
//  bool have2Diso;
////  #endif //GEOMETRY == CARTESIAN or GEOMETRY == SPHERICAL
////  #if GEOMETRY == SPHERICAL
//  bool haveVsh;
  ForcingType forcingType;
  int normal2Diso;

  IdefixArray2D<real> k3Diso;
  IdefixHostArray2D<real> k3DisoHost;
  real kmin;
  real kmax;
  int ellmin;
  int ellmax;
  int mmin;
  int mmax;
//  #endif //GEOMETRY == SPHERICAL
};

#endif // FORCING_FORCING_HPP_
