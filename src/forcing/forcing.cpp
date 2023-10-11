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
    this->lmin = input.GetOrSet<int>("Forcing","lmin",0, 0);
    this->mmin = input.GetOrSet<int>("Forcing","mmin",0, 0);
    this->lmax = input.GetOrSet<int>("Forcing","lmax",0, 3);
    this->mmax = input.GetOrSet<int>("Forcing","mmax",0, 3);
    if (mmax > lmax) {
      IDEFIX_ERROR("mmax cannot be greater than lmax");
    }
    if (mmin > lmin) {
      IDEFIX_ERROR("mmin cannot be greater than lmin");
    }
    this->t_corr = input.GetOrSet<real>("Forcing","t_corr",0, 1.);
    this->eps_Ylm = input.GetOrSet<real>("Forcing","eps_Ylm",0, 1.);
    this->eps_Slm = input.GetOrSet<real>("Forcing","eps_Slm",0, 1.);
    this->eps_Tlm = input.GetOrSet<real>("Forcing","eps_Tlm",0, 1.);

    // Allocate required arrays
    this->forcingTerm = IdefixArray4D<real>("ForcingTerm", COMPONENTS,
                                data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);
  
    this->OUprocesses.InitProcesses(this->lmin,this->lmax,this->mmin,this->mmax,0.,t_corr,eps_Ylm,eps_Slm,eps_Tlm);
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
    idfx::cout << "Forcing: lmin=" << lmin << " and lmax=" << lmax << "." << std::endl;
    idfx::cout << "Forcing: mmin=" << mmin << " and mmax=" << mmax << "." << std::endl;
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
    int lmin = this->lmin;
    int lmax = this->lmax;
    int mmin = this->mmin;
    int mmax = this->mmax;
    IdefixArray2D<real> jl = this->data->jl;
    IdefixArray4D<real> Ylm_r = this->data->Ylm_r;
    IdefixArray4D<real> Slm_th = this->data->Slm_th;
    IdefixArray4D<real> Slm_phi = this->data->Slm_phi;
    IdefixArray4D<real> Tlm_th = this->data->Tlm_th;
    IdefixArray4D<real> Tlm_phi = this->data->Tlm_phi;
    IdefixArray4D<real> forcingTerm = this->forcingTerm;
    IdefixArray3D<real> OuValues = this->OUprocesses.ouValues;
    idefix_for("ComputeForcing", 0, data->np_tot[KDIR], 0, data->np_tot[JDIR], 0, data->np_tot[IDIR],
                KOKKOS_LAMBDA (int k, int j, int i) {
                    real forcing_MX1=0;
                    real forcing_MX2=0;
                    real forcing_MX3=0;
                    for (int l=lmin; l<lmax; l++) {
//                      real jli = jl(l,i);
                      real jli = 1.;
                      for (int m=mmin; m<mmax; m++) {
                        forcing_MX1 += jli*Ylm_r(l,m,k,j)*OuValues(0,l,m);
                        forcing_MX2 += jli*Slm_th(l,m,k,j)*OuValues(1,l,m) + jli*Tlm_th(l,m,k,j)*OuValues(2,l,m);
                        forcing_MX3 += jli*Slm_phi(l,m,k,j)*OuValues(1,l,m) + jli*Tlm_phi(l,m,k,j)*OuValues(2,l,m);
                      }
                    }
   
                    forcingTerm(IDIR,k,j,i) += forcing_MX1;
                    forcingTerm(JDIR,k,j,i) += forcing_MX2;
                    forcingTerm(KDIR,k,j,i) += forcing_MX3;
    });
  #endif // GEOMETRY == SPHERICAL

  OUprocesses.UpdateProcessesValues(dt);
//  OUprocesses.UpdateProcessesValues(dt, eps_Ylm, eps_Slm, eps_Tlm);
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

