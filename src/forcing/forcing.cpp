// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#include <string>

#include "forcing.hpp"
#include "dataBlock.hpp"
#include "dataBlockHost.hpp"
#include "input.hpp"

Forcing::Forcing(Input &input, DataBlock *datain) {
  idfx::pushRegion("Forcing::Forcing");
  this->data = datain;

  // Forcing
  #if GEOMETRY == SPHERICAL
    if(input.CheckEntry("Forcing","lmax")>=0) {
      this->lmax = input.Get<int>("Forcing","lmax",0);
    } else {
      idfx::cout << "No lmax provided. Assuming 2." << std::endl;
      this->lmax = 2;
    }
    if(input.CheckEntry("Forcing","mmax")>=0) {
      this->mmax = input.Get<int>("Forcing","mmax",0);
    } else {
      idfx::cout << "No mmax provided. Assuming 2." << std::endl;
      this->mmax = 2;
    }
    if(input.CheckEntry("Forcing","t_corr")>=0) {
      this->t_corr = input.Get<real>("Forcing","t_corr",0);
    } else {
      idfx::cout << "No correlation time provided. Assuming 1." << std::endl;
      this->t_corr = 1.;
    }
    if(input.CheckEntry("Forcing","frms_Ylm")>=0) {
      this->frms_Ylm = input.Get<real>("Forcing","frms_Ylm",0);
    } else {
      idfx::cout << "No rms forcing for Ylm provided. Assuming 1." << std::endl;
      this->frms_Ylm = 1.;
    }
    if(input.CheckEntry("Forcing","frms_Slm")>=0) {
      this->frms_Slm = input.Get<real>("Forcing","frms_Slm",0);
    } else {
      idfx::cout << "No rms forcing for Slm provided. Assuming 1." << std::endl;
      this->frms_Slm = 1.;
    }
    if(input.CheckEntry("Forcing","frms_Tlm")>=0) {
      this->frms_Tlm = input.Get<real>("Forcing","frms_Tlm",0);
    } else {
      idfx::cout << "No rms forcing for Tlm provided. Assuming 1." << std::endl;
      this->frms_Tlm = 1.;
    }
  #endif // GEOMETRY == SPHERICAL

  // Allocate required arrays
  this->forcingTerm = IdefixArray4D<real>("ForcingTerm", COMPONENTS,
                              data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);

  int ntheta = input.Get<int>("Grid","X2-grid",2); // number of points in the latitude direction  (constraint: nlat >= lmax+1)
  int nphi = input.Get<int>("Grid","X3-grid",2);  // number of points in the longitude direction (constraint: nphi >= 2*mmax+1)
  int nr_proc = data->np_int[IDIR];
  int ntheta_proc = data->np_int[JDIR];
  int nphi_proc = data->np_int[KDIR];
  int koffset = data->gbeg[KDIR] - data->nghost[KDIR];
  int joffset = data->gbeg[JDIR] - data->nghost[JDIR];
  DataBlockHost d(*datain, 1);
  this->vsh = std::make_unique<Vsh>(&d);
  this->vsh->GenerateCellVsh(0);
  this->vsh->Generatejl();

  this->OUprocesses.InitProcesses(lmax,mmax,0.,t_corr);

//  this->skipGravity = input.GetOrSet<int>("Gravity","skip",0,1);
//  if(skipGravity<1) {
//    IDEFIX_ERROR("[Gravity]:skip should be a strictly positive integer");
//  }
  idfx::popRegion();
}

void Forcing::ShowConfig() {
  #if GEOMETRY == SPHERICAL
    idfx::cout << "Forcing: ENABLED." << std::endl;
    idfx::cout << "Forcing: lmax=" << lmax << " and mmax=" << mmax << " ." << std::endl;
    idfx::cout << "Forcing: t_corr=" << this->t_corr << " ." << std::endl;
    idfx::cout << "Forcing: frms_Ylm=" << frms_Ylm << ", frms_Slm=" << frms_Slm << ", frms_Tlm =" << frms_Tlm << " ." << std::endl;
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
                      real jli = vsh->jl(l,i);
                      for (int m=0; m<this->mmax; m++) {
                        forcing_MX1 += jli*this->vsh->Ylm_r(l,m,k,j)*this->OUprocesses.ouValues(0,l,m);
                        forcing_MX2 += jli*this->vsh->Slm_th(l,m,k,j)*this->OUprocesses.ouValues(1,l,m) + jli*this->vsh->Tlm_th(l,m,k,j)*this->OUprocesses.ouValues(2,l,m);
                        forcing_MX3 += jli*vsh->Slm_phi(l,m,k,j)*this->OUprocesses.ouValues(1,l,m) + jli*this->vsh->Tlm_phi(l,m,k,j)*this->OUprocesses.ouValues(2,l,m);
                      }
                    }
   
                    this->forcingTerm(IDIR,k,j,i) += forcing_MX1;
                    this->forcingTerm(JDIR,k,j,i) += forcing_MX2;
                    this->forcingTerm(KDIR,k,j,i) += forcing_MX3;
      });
    });
  #endif // GEOMETRY == SPHERICAL

  OUprocesses.UpdateProcesses(dt, frms_Ylm, frms_Slm, frms_Tlm);
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

