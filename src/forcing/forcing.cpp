// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#include <string>

#include "forcing.hpp"
#include "dataBlock.hpp"
#include "input.hpp"

#if GEOMETRY == SPHERICAL
  #include "vsh.hpp"
#endif

Forcing::Forcing(Input &input, DataBlock *datain) {
  idfx::pushRegion("Forcing::Forcing");
  this->data = datain;

  // Forcing
  if(input.CheckEntry("Forcing","OU")>=0) {
    std::string forcingString = input.Get<std::string>("Forcing","OU",0);
    // GET BACK INFO
  }

  // Allocate required arrays
  this->forcingTerm = IdefixArray4D<real>("ForcingTerm", COMPONENTS,
                              data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);

//  this->skipGravity = input.GetOrSet<int>("Gravity","skip",0,1);
//  if(skipGravity<1) {
//    IDEFIX_ERROR("[Gravity]:skip should be a strictly positive integer");
//  }
  idfx::popRegion();
}

void Forcing::ShowConfig() {
  idfx::cout << "Forcing: ENABLED with t_corr=." << this->t_corr << " ." << std::endl;
//    if(skipGravity>1) {
//      idfx::cout << "Gravity: gravity field will be updated every " << skipGravity
//                 << " cycles." << std::endl;
//    }
  
}
// This function compute the gravitational field, using both body force and potential
void Forcing::ComputeForcing(int stepNumber) {
  idfx::pushRegion("Forcing::ComputeForcing");
  // DO STUFF
  idfx::popRegion();
}

//// Fill the gravitational potential with zeros
//void Gravity::ResetPotential() {
//  idfx::pushRegion("Gravity::ResetPotential");
//  IdefixArray3D<real> phiP = this->phiP;
//  idefix_for("Gravity::ResetPotential",
//              0, data->np_tot[KDIR],
//              0, data->np_tot[JDIR],
//              0, data->np_tot[IDIR],
//              KOKKOS_LAMBDA(int k, int j, int i) {
//                phiP(k,j,i) = ZERO_F;
//              });
//  idfx::popRegion();
//}

