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

Forcing::Forcing(Input &input, DataBlock *datain) {
  idfx::pushRegion("Forcing::Forcing");
  this->data = datain;

  // Forcing
  #if GEOMETRY == SPHERICAL
    this->lmax = input.GetOrSet<int>("Forcing","lmax",0, 3);
    this->mmax = input.GetOrSet<int>("Forcing","mmax",0, 3);
    this->t_corr = input.GetOrSet<real>("Forcing","t_corr",0, 1.);
    this->eps_Ylm = input.GetOrSet<real>("Forcing","eps_Ylm",0, 1.);
    this->eps_Slm = input.GetOrSet<real>("Forcing","eps_Slm",0, 1.);
    this->eps_Tlm = input.GetOrSet<real>("Forcing","eps_Tlm",0, 1.);

    // Allocate required arrays
    this->forcingTerm = IdefixArray4D<real>("ForcingTerm", COMPONENTS,
                                data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);
  
    this->OUprocesses.InitProcesses(this->lmax,this->mmax,0.,t_corr);
  #endif // GEOMETRY == SPHERICAL

//  this->skipGravity = input.GetOrSet<int>("Gravity","skip",0,1);
//  if(skipGravity<1) {
//    IDEFIX_ERROR("[Gravity]:skip should be a strictly positive integer");
//  }
  idfx::popRegion();
}

void Forcing::ShowConfig() {
  #if GEOMETRY == SPHERICAL
    idfx::cout << "Forcing: ENABLED." << std::endl;
    idfx::cout << "Forcing: lmax=" << lmax << " and mmax=" << mmax << "." << std::endl;
    idfx::cout << "Forcing: t_corr=" << this->t_corr << " ." << std::endl;
    idfx::cout << "Forcing: eps_Ylm=" << eps_Ylm << ", eps_Slm=" << eps_Slm << ", eps_Tlm =" << eps_Tlm << " ." << std::endl;
  #endif // GEOMETRY == SPHERICAL
//    if(skipGravity>1) {
//      idfx::cout << "Gravity: gravity field will be updated every " << skipGravity
//                 << " cycles." << std::endl;
//    }
  
}
// This function compute the acceleration part of the forcing term
//void Forcing::ComputeForcing(int stepNumber) {
void Forcing::ComputeForcing(real dt) {
  idfx::pushRegion("Forcing::ComputeForcing");
  
  ResetForcingTerm();
  IdefixArray1D<real> x1 = data->x[IDIR];

  #if GEOMETRY == CARTESIAN
    IDEFIX_ERROR("Forcing in CARTESIAN geometry not implemented");
  #endif // GEOMETRY == CARTESIAN
  #if GEOMETRY == POLAR
    IDEFIX_ERROR("Forcing in POLAR geometry not implemented");
  #endif // GEOMETRY == POLAR
  #if GEOMETRY == CYLINDRICAL
    IDEFIX_ERROR("Forcing in CYLINDRICAL geometry not implemented");
  #endif // GEOMETRY == CYLINDRICAL
  #if GEOMETRY == SPHERICAL
//    idefix_for("ComputeForcing", 0, data->np_tot[KDIR], 0, data->np_tot[JDIR], 0, data->np_tot[IDIR],
//                KOKKOS_LAMBDA (int k, int j, int i) {
    idefix_for("ComputeForcingThetaPhi", 0, data->np_tot[KDIR], 0, data->np_tot[JDIR],
                KOKKOS_LAMBDA (int k, int j) {
      idefix_for("ComputeForcingR", 0, data->np_tot[IDIR],
                  KOKKOS_LAMBDA (int i) {
                    real forcing_MX1=0;
                    real forcing_MX2=0;
                    real forcing_MX3=0;
                    for (int l=0; l<this->lmax; l++) {
                      real jli = this->data->jl(l,i);
                      for (int m=0; m<this->mmax; m++) {
                        forcing_MX1 += jli*this->data->Ylm_r(l,m,k,j)*this->OUprocesses.ouValues(0,l,m);
                        forcing_MX2 += jli*this->data->Slm_th(l,m,k,j)*this->OUprocesses.ouValues(1,l,m) + jli*this->data->Tlm_th(l,m,k,j)*this->OUprocesses.ouValues(2,l,m);
                        forcing_MX3 += jli*data->Slm_phi(l,m,k,j)*this->OUprocesses.ouValues(1,l,m) + jli*this->data->Tlm_phi(l,m,k,j)*this->OUprocesses.ouValues(2,l,m);
                      }
                    }
   
                    this->forcingTerm(IDIR,k,j,i) += forcing_MX1;
                    this->forcingTerm(JDIR,k,j,i) += forcing_MX2;
                    this->forcingTerm(KDIR,k,j,i) += forcing_MX3;
      });
    });
  #endif // GEOMETRY == SPHERICAL

  OUprocesses.UpdateProcesses(dt, eps_Ylm, eps_Slm, eps_Tlm);
  idfx::popRegion();
}

// Fill the forcing term with zeros
void Forcing::ResetForcingTerm() {
  idfx::pushRegion("Forcing::ResetForcingTerm");
  IdefixArray4D<real> forcingTerm = this->forcingTerm;
  idefix_for("Forcing::ResetForcingTerm",
              0, data->np_tot[KDIR],
              0, data->np_tot[JDIR],
              0, data->np_tot[IDIR],
              KOKKOS_LAMBDA(int k, int j, int i) {
                forcingTerm(IDIR,k,j,i) = ZERO_F;
                forcingTerm(JDIR,k,j,i) = ZERO_F;
                forcingTerm(KDIR,k,j,i) = ZERO_F;
              });
  idfx::popRegion();
}

