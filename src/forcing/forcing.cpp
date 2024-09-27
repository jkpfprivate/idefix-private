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
  this->seed = input.GetOrSet<int>("Forcing","seed",0,0);

  // Forcing
//  #if GEOMETRY == SPHERICAL & VSH == YES
//    this->lmin = input.GetOrSet<int>("Forcing","lmin",0, 0);
//    this->mmin = input.GetOrSet<int>("Forcing","mmin",0, 0);
//    this->lmax = input.GetOrSet<int>("Forcing","lmax",0, 3);
//    this->mmax = input.GetOrSet<int>("Forcing","mmax",0, 3);
//    if (mmax > lmax) {
//      IDEFIX_ERROR("Forcing: mmax cannot be greater than lmax");
//    } else if (mmin > lmin) {
//      IDEFIX_ERROR("Forcing: mmin cannot be greater than lmin");
//    } else if (lmin >= lmax) {
//      IDEFIX_ERROR("Forcing: lmax cannot be less than or equal to lmax");
//    } else if (mmin >= mmax) {
//      IDEFIX_ERROR("Forcing: mmin cannot be less than or equal to mmax");
//    }
//    this->t_corr = input.GetOrSet<real>("Forcing","t_corr",0, 1.);
//    this->eps_Ylm = input.GetOrSet<real>("Forcing","eps_Ylm",0, 1.);
//    this->eps_Slm = input.GetOrSet<real>("Forcing","eps_Slm",0, 1.);
//    this->eps_Tlm = input.GetOrSet<real>("Forcing","eps_Tlm",0, 1.);

  this->write = input.GetOrSet<int>("Forcing","write",0, 0);
  std::string folder = input.GetOrSet<std::string>("Forcing","filename",0, "testOU");

//  #endif // GEOMETRY == SPHERICAL & VSH = YES

//  this->have3Diso = false;
//  this->have2Diso = false;
//  this->haveVsh = false;
//  this->haveUserDef = false;

//  this->nForcingModesIdir = 0;
//  this->nForcingModesJdir = 0;
//  this->nForcingModesKdir = 0;
  this->nForcingModes = 0;

  this->kmin = -1.;
  this->kmax = -1.;
  this->ellmin = -1.;
  this->ellmax = -1.;
  this->mmin = -1.;
  this->mmax = -1.;

  real kx0, ky0, kz0;

  if (input.CheckEntry("Forcing", "3Diso")>=0) {
//    this->have3Diso = true;
    this->forcingType = iso3D;
    #if GEOMETRY != CARTESIAN
      IDEFIX_ERROR("You cannot have 3D isotropic forcing in another geometry than Cartesian");
    #endif //GEOMETRY != CARTESIAN
    #if COMPONENTS < 3 or DIMENSIONS < 3
      IDEFIX_ERROR("You cannot have 3D isotropic forcing with less than 3 components and dimensions.");
    #endif //COMPONENTS < 3 or DIMENSIONS < 3
    this->kmin = input.Get<int>("Forcing","3Diso", 1);
    this->kmax = input.Get<int>("Forcing","3Diso", 2);

//    k3Diso = std::vector<std::vector<real>>(0);
    std::vector<std::vector<real>> k3Disovec;

    kx0 = 2.*M_PI/(data->mygrid->xend[IDIR] - data->mygrid->xbeg[IDIR]);
    ky0 = 2.*M_PI/(data->mygrid->xend[JDIR] - data->mygrid->xbeg[JDIR]);
    kz0 = 2.*M_PI/(data->mygrid->xend[KDIR] - data->mygrid->xbeg[KDIR]);
    int nxmax = kmax/kx0;
    int nymax = kmax/ky0;
    int nzmax = kmax/kz0;
//    int nxmin = kmin/kx0;
//    int nymin = kmin/ky0;
//    int nzmin = kmin/kz0;
//    for (int nx=nxmin; nx<nxmax; nx++) {
//      for (int ny=nymin; ny<nymax; ny++) {
//        for (int nz=nzmin; nz<nzmax; nz++) {
    for (int nx=0; nx<nxmax; nx++) {
      for (int ny=0; ny<nymax; ny++) {
        for (int nz=0; nz<nzmax; nz++) {
          real kx = kx0*nx;
          real ky = ky0*ny;
          real kz = kz0*nz;
          real k_2 = kx*kx + ky*ky + kz*kz;
          if (k_2 >= kmin*kmin and k_2 <= kmax*kmax) {
            nForcingModes ++;
            k3Disovec.push_back({kx, ky, kz});
          }
        }
      }
    }
    k3DisoHost = IdefixHostArray2D<real>("k3DisoHost", nForcingModes);
    k3Diso = IdefixArray2D<real>("k3Diso", nForcingModes, 3);
    for (int l=0; l<nForcingModes; l++) {
      k3DisoHost(l,0) = k3Disovec[l][0];
      k3DisoHost(l,1) = k3Disovec[l][1];
      k3DisoHost(l,2) = k3Disovec[l][2];
    }
    Kokkos::deep_copy(k3Diso, k3DisoHost);
   

    this->forcingModes = IdefixArray4D<Kokkos::complex<real>>("forcingModes", nForcingModes);
    this->forcingModesIdir = IdefixArray4D<real>("forcingModesIdir", nForcingModes);
    this->forcingModesJdir = IdefixArray4D<real>("forcingModesJdir", nForcingModes);
    this->forcingModesKdir = IdefixArray4D<real>("forcingModesKdir", nForcingModes);
  }

  else if (input.CheckEntry("Forcing", "2Diso")>=0) {
//    this->have2Diso = true;
    this->forcingType = iso2D;
    this->normal2Diso = input.Get<int>("Forcing","3Diso", 1);
    if (normal2Diso < 0 or normal2Diso > 2) IDEFIX_ERROR("The normal component for 2D isotropic forcing cannot not be something else than IDIR, JDIR or KDIR");
    #if COMPONENTS < 2 or DIMENSIONS < 2
      IDEFIX_ERROR("You cannot have 2D isotropic forcing with less than 2 components and dimensions.");
    #endif //COMPONENTS < 3 or DIMENSIONS < 3
    this->kmin = input.Get<int>("Forcing","3Diso", 2);
    this->kmax = input.Get<int>("Forcing","3Diso", 3);
  }

  else if (input.CheckEntry("Forcing", "vsh")>=0) {
//    this->haveVsh = true;
    this->forcingType = vsh;
    #if GEOMETRY != SPHERICAL
      IDEFIX_ERROR("You cannot have vsh forcing in another geometry than spherical.");
    #endif //GEOMETRY != SPHERICAL
    #if VSH == NO
      IDEFIX_ERROR("You cannot have vsh forcing without the VSH module.");
    #endif //VSH == NO
    this->ellmin = input.Get<int>("Forcing","vsh", 1);
    this->ellmax = input.Get<int>("Forcing","vsh", 2);
    this->mmin = input.GetOrSet<int>("Forcing","vsh", 3, 0);
    this->mmax = input.GetOrSet<int>("Forcing","vsh", 4, this->ellmax);
    for (int ell=ellmin; ell<ellmax; ell++) {
      for (int m=mmin; m<mmax & m<ell+1 ; m++) {
        this->nForcingModesIdir++;
        this->nForcingModesJdir++;
        this->nForcingModesKdir++;
      }
    }
    // WARNING, vsh case needs a special treatment since the forcing modes between
    // the theta and phi directions are coupled to each other...
    this->nForcingModesJdir *= 2;
    this->nForcingModesKdir *= 2;
    #if VSH == NO
      IDEFIX_ERROR("You cannot have vsh forcing without the VSH module.");
    #endif //VSH == NO
    this->forcingModesIdir = IdefixArray4D<real>("forcingModesIdir", nForcingModesIdir);
    this->forcingModesJdir = IdefixArray4D<real>("forcingModesJdir", nForcingModesJdir);
    this->forcingModesKdir = IdefixArray4D<real>("forcingModesKdir", nForcingModesKdir);
  }

//  else { this->haveUserDef = true;}
  else { this->forcingType = userDef;}

  this->nForcingModes = this->nForcingModesIdir + this->nForcingModesJdir + this->nForcingModesKdir;

  // Allocate required arrays
  this->forcingTerm = IdefixArray4D<real>("ForcingTerm", COMPONENTS,
                              data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);
  this->means = IdefixArray1D<real>("means", nForcingModes);
  this->tcorrs = IdefixArray1D<real>("tcorrs", nForcingModes);
  this->epsilons = IdefixArray1D<real>("epsilons", nForcingModes);

  this->OUprocesses.InitProcesses(folder, this->seed, this->nForcingModes, this->means, this->tcorrs, this->epsilons);

  idfx::popRegion();
}

void Forcing::ShowConfig() {
  idfx::cout << "Forcing: ENABLED with seed " << seed << "." << std::endl;
  switch(forcingType) {
    case ForcingType::iso3D:
      idfx::cout << "Forcing: 3D isotropic." << std::endl;
      idfx::cout << "Forcing: kmin=" << kmin << " and kmax" << kmax << " ." << std::endl;
      break;
    case ForcingType::iso2D:
      idfx::cout << "Forcing: 2D isotropic with normal." << normal2Diso << "." << std::endl;
      idfx::cout << "Forcing: kmin=" << kmin << " and kmax" << kmax << " ." << std::endl;
      break;
    case ForcingType::vsh:
      idfx::cout << "Forcing: vector spherical harmonics." << std::endl;
      idfx::cout << "Forcing: ellmin=" << ellmin << " and ellmax" << ellmax << " ." << std::endl;
      break;
    case ForcingType::userDef:
      idfx::cout << "Forcing: userdef." << std::endl;
      break;
  }

//    if(skipGravity>1) {
//      idfx::cout << "Gravity: gravity field will be updated every " << skipGravity
//                 << " cycles." << std::endl;
//    }
}

void Forcing::InitForcingModes() {
  idfx::pushRegion("Forcing::InitForcingModes");
  IdefixArray4D<Kokkos::complex<real>> forcingModes = this->forcingModes;
  IdefixArray4D<real> forcingModesIdir = this->forcingModesIdir;
  IdefixArray4D<real> forcingModesJdir = this->forcingModesJdir;
  IdefixArray4D<real> forcingModesKdir = this->forcingModesKdir;
  IdefixArray1D<real> x1 = data->x[IDIR];
  IdefixArray1D<real> x2 = data->x[JDIR];
  IdefixArray1D<real> x3 = data->x[KDIR];
  IdefixArray2D<real> k3Diso = this->k3Diso;
  Kokkos::complex unit_j(0.,1.);
  idefix_for("Forcing::InitForcingModes",
              0, nForcingModes,
              0, data->np_tot[KDIR],
              0, data->np_tot[JDIR],
              0, data->np_tot[IDIR],
              KOKKOS_LAMBDA(int l, int k, int j, int i) {
                forcingModesIdir(l,k,j,i) = ZERO_F;
                forcingModesJdir(l,k,j,i) = ZERO_F;
                forcingModesKdir(l,k,j,i) = ZERO_F;
              });
  real kx0, ky0, kz0;
  switch(forcingType) {
    case ForcingType::iso3D:
      kx0 = 2.*M_PI/(data->mygrid->xend[IDIR] - data->mygrid->xbeg[IDIR]);
      ky0 = 2.*M_PI/(data->mygrid->xend[JDIR] - data->mygrid->xbeg[JDIR]);
      kz0 = 2.*M_PI/(data->mygrid->xend[KDIR] - data->mygrid->xbeg[KDIR]);
      idefix_for("ComputeForcing", 0, nForcingModes, 0, data->np_tot[KDIR], 0, data->np_tot[JDIR], 0, data->np_tot[IDIR],
                  KOKKOS_LAMBDA (int l, int k, int j, int i) {
                    real kx = k3Diso(l, IDIR);
                    real ky = k3Diso(l, JDIR);
                    real kz = k3Diso(l, KDIR);
//                    real kx = k3Diso[l][IDIR];
//                    real ky = k3Diso[l][JDIR];
//                    real kz = k3Diso[l][KDIR];
                    forcingModes(l,k,j,i) = exp(unit_j*(kx*x1(i) + ky*x2(j) + ky*x3(k)));
      });
      break;
    case ForcingType::iso2D:
      kx0 = 2.*M_PI/(data->mygrid->xend[IDIR] - data->mygrid->xbeg[IDIR]);
      ky0 = 2.*M_PI/(data->mygrid->xend[JDIR] - data->mygrid->xbeg[JDIR]);
      kz0 = 2.*M_PI/(data->mygrid->xend[KDIR] - data->mygrid->xbeg[KDIR]);
      break;
    #if VSH == YES
      case ForcingType::vsh:
        int ellmin = this->ellmin;
        int ellmax = this->ellmax;
        int mmin = this->mmin;
        int mmax = this->mmax;
        IdefixArray4D<real> Ylm_r = this->data->Ylm_r;
        IdefixArray4D<real> Slm_th = this->data->Slm_th;
        IdefixArray4D<real> Slm_phi = this->data->Slm_phi;
        IdefixArray4D<real> Tlm_th = this->data->Tlm_th;
        IdefixArray4D<real> Tlm_phi = this->data->Tlm_phi;
        idefix_for("ComputeForcing", 0, data->np_tot[KDIR], 0, data->np_tot[JDIR], 0, data->np_tot[IDIR],
                    KOKKOS_LAMBDA (int k, int j, int i) {
                      for (int ell=ellmin; ell<ellmax; ell++) {
                        for (int m=mmin; m<mmax & m<ell+1 ; m++) {
                          forcingModesIdir(ell,k,j,i) = Ylm_r(ell,m,j,i);
                        }
                      }
        });
        break;
    #endif //VSH == YES
//    case ForcingType::userDef:
//      break;
  }
}

// This function compute the acceleration part of the forcing term
//void Forcing::ComputeForcing(int stepNumber) {
void Forcing::ComputeForcing(real dt) {
  idfx::pushRegion("Forcing::ComputeForcing");
  
  ResetForcingTerm();

//  #if GEOMETRY == CARTESIAN
//    IDEFIX_ERROR("Forcing in CARTESIAN geometry not implemented");
//  #endif // GEOMETRY == CARTESIAN
//  #if GEOMETRY == POLAR
//    IDEFIX_ERROR("Forcing in POLAR geometry not implemented");
//  #endif // GEOMETRY == POLAR
//  #if GEOMETRY == CYLINDRICAL
//    IDEFIX_ERROR("Forcing in CYLINDRICAL geometry not implemented");
//  #endif // GEOMETRY == CYLINDRICAL
//  #if GEOMETRY == SPHERICAL
//    #if VSH == YES
  IdefixArray4D<real> forcingModesIdir = this->forcingModesIdir;
  IdefixArray4D<real> forcingModesJdir = this->forcingModesJdir;
  IdefixArray4D<real> forcingModesKdir = this->forcingModesKdir;
  IdefixArray4D<real> forcingTerm = this->forcingTerm;
  IdefixArray1D<real> OuValues = this->OUprocesses.ouValues;
  int nForcingModes = this->nForcingModes;
  idefix_for("ComputeForcing", 0, data->np_tot[KDIR], 0, data->np_tot[JDIR], 0, data->np_tot[IDIR],
              KOKKOS_LAMBDA (int k, int j, int i) {
                  real forcing_MX1=0;
                  real forcing_MX2=0;
                  real forcing_MX3=0;
                  for (int l=0; l<nForcingModes; l++) {
                    forcing_MX1 += OuValues(l)*forcingModesIdir(l,k,j,i);
                    #if COMPONENTS >= 2
                      forcing_MX2 += OuValues(l)*forcingModesJdir(l,k,j,i);
                      #if COMPONENTS == 3
                        forcing_MX3 += OuValues(l)*forcingModesKdir(l,k,j,i);
                      #endif //COMPONENTS == 3
                    #endif //COMPONENTS >= 2
                  }
                  forcingTerm(IDIR,k,j,i) += real(forcing_MX1);
                  forcingTerm(JDIR,k,j,i) += real(forcing_MX2);
                  forcingTerm(KDIR,k,j,i) += real(forcing_MX3);
  });
//    #endif // VSH == YES
//  #endif // GEOMETRY == SPHERICAL

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

