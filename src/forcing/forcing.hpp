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
  IdefixArray4D<Kokkos::complex<real>> forcingModesIdir;
  #if COMPONENTS >= 2
    IdefixArray4D<Kokkos::complex<real>> forcingModesJdir;
    #if COMPONENTS == 3
      IdefixArray4D<Kokkos::complex<real>> forcingModesKdir;
    #endif //COMPONENTS == 3
  #endif //COMPONENTS >= 2
  // Forcing term
  IdefixArray4D<real> forcingTerm;
  OrnsteinUhlenbeckProcesses oUprocesses;

//  // Whether we should skip gravity computation every n steps
//  int skipGravity{1};

  int write;
 private:

  DataBlock *data;
  int seed;

  int nForcingModes;
  IdefixArray2D<real> tcorrs;
  IdefixArray2D<real> means;
  IdefixArray2D<real> epsilons;

  ForcingType forcingType;
  int normal2Diso;

  IdefixArray2D<real> k3Diso;
  IdefixHostArray2D<real> k3DisoHost;
  IdefixArray2D<real> k2Diso;
  IdefixHostArray2D<real> k2DisoHost;
  IdefixArray2D<real> ellmVsh;
  IdefixHostArray2D<real> ellmVshHost;
  real kmin;
  real kmax;
  int ellmin;
  int ellmax;
  int mmin;
  int mmax;
//  #endif //GEOMETRY == SPHERICAL
};

#endif // FORCING_FORCING_HPP_
