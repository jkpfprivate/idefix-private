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

  this->write = input.GetOrSet<int>("Forcing","write",0, 0);
  std::string folder = input.GetOrSet<std::string>("Forcing","filename",0, "testOU");

  this->nForcingModes = 0;

  this->kmin = -1.;
  this->kmax = -1.;
  this->ellmin = -1.;
  this->ellmax = -1.;
  this->mmin = -1.;
  this->mmax = -1.;
  this->haveSolenoidalForcing = false;

  real kx0, ky0, kz0;

  if (input.CheckEntry("Forcing", "iso3D")>=0) {
    this->forcingType = iso3D;
    #if GEOMETRY != CARTESIAN
      IDEFIX_ERROR("You cannot have 3D isotropic forcing in another geometry than Cartesian");
    #endif //GEOMETRY != CARTESIAN
    #if COMPONENTS < 3 or DIMENSIONS < 3
      IDEFIX_ERROR("You cannot have 3D isotropic forcing with less than 3 components and dimensions.");
    #endif //COMPONENTS < 3 or DIMENSIONS < 3
// if () IDEFIX_ERROR("You should not have 3D isotropic forcing with other boundary conditions than fully periodic.");
    #if GEOMETRY == SPHERICAL
      IDEFIX_ERROR("Comparing radial with angular wavenumbers is somewhat shady.");
    #endif //GEOMETRY == SPHERICAL
    this->kmin = input.Get<real>("Forcing","iso3D", 0);
    this->kmax = input.Get<real>("Forcing","iso3D", 1);

    std::vector<std::vector<real>> k3Disovec;

    kx0 = 2.*M_PI/(data->mygrid->xend[IDIR] - data->mygrid->xbeg[IDIR]);
    ky0 = 2.*M_PI/(data->mygrid->xend[JDIR] - data->mygrid->xbeg[JDIR]);
    kz0 = 2.*M_PI/(data->mygrid->xend[KDIR] - data->mygrid->xbeg[KDIR]);
    int nxmax = kmax/kx0 + 1;
    int nymax = kmax/ky0 + 1;
    int nzmax = kmax/kz0 + 1;
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
            modeNames.push_back({"I" + std::to_string(nx) + std::to_string(ny) + std::to_string(nz), "J" + std::to_string(nx) + std::to_string(ny) + std::to_string(nz), "K" + std::to_string(nx) + std::to_string(ny) + std::to_string(nz)});
          }
        }
      }
    }
    k3DisoHost = IdefixHostArray2D<real>("k3DisoHost", nForcingModes, 3);
    k3Diso = IdefixArray2D<real>("k3Diso", nForcingModes, 3);
    for (int l=0; l<nForcingModes; l++) {
      k3DisoHost(l,0) = k3Disovec[l][0];
      k3DisoHost(l,1) = k3Disovec[l][1];
      k3DisoHost(l,2) = k3Disovec[l][2];
    }
    Kokkos::deep_copy(k3Diso, k3DisoHost);
    this->forcingModesIdir = IdefixArray4D<Kokkos::complex<real>>("forcingModesIdir", nForcingModes,
                              data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);
    this->forcingModesJdir = IdefixArray4D<Kokkos::complex<real>>("forcingModesJdir", nForcingModes,
                              data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);
    this->forcingModesKdir = IdefixArray4D<Kokkos::complex<real>>("forcingModesKdir", nForcingModes,
                              data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);
  }

  else if (input.CheckEntry("Forcing", "iso2D")>=0) {
    this->forcingType = iso2D;
    this->normal2DisoStr = input.Get<std::string>("Forcing","iso2D", 0);
    if (this->normal2DisoStr == "IDIR") this->normal2Diso = IDIR;
    else if (this->normal2DisoStr == "JDIR") this->normal2Diso = JDIR;
    else if (this->normal2DisoStr == "KDIR") this->normal2Diso = KDIR;
    else IDEFIX_ERROR("The normal component for 2D isotropic forcing cannot not be something else than IDIR, JDIR or KDIR");
    #if COMPONENTS < 2 or DIMENSIONS < 2
      IDEFIX_ERROR("You cannot have 2D isotropic forcing with less than 2 components and dimensions.");
    #endif //COMPONENTS < 2 or DIMENSIONS < 2
    #if GEOMETRY == SPHERICAL
      if (this->normal2DisoStr == "JDIR" or this->normal2DisoStr == "KDIR") IDEFIX_ERROR("Comparing radial with angular wavenumbers is somewhat shady.");
    #endif //GEOMETRY == SPHERICAL

    this->kmin = input.Get<real>("Forcing","iso2D", 1);
    this->kmax = input.Get<real>("Forcing","iso2D", 2);

    std::vector<std::vector<real>> k2Disovec;

    kx0 = (this->normal2Diso==IDIR) ? ZERO_F : 2.*M_PI/(data->mygrid->xend[IDIR] - data->mygrid->xbeg[IDIR]);
    ky0 = (this->normal2Diso==JDIR) ? ZERO_F : 2.*M_PI/(data->mygrid->xend[JDIR] - data->mygrid->xbeg[JDIR]);
    kz0 = (this->normal2Diso==KDIR) ? ZERO_F : 2.*M_PI/(data->mygrid->xend[KDIR] - data->mygrid->xbeg[KDIR]);
    int nxmax = (this->normal2Diso==IDIR) ? 1 : kmax/kx0 + 1;
    int nymax = (this->normal2Diso==JDIR) ? 1 : kmax/ky0 + 1;
    int nzmax = (this->normal2Diso==KDIR) ? 1 : kmax/kz0 + 1;
    for (int nx=0; nx<nxmax; nx++) {
      for (int ny=0; ny<nymax; ny++) {
        for (int nz=0; nz<nzmax; nz++) {
          real kx = kx0*nx;
          real ky = ky0*ny;
          real kz = kz0*nz;
          real k_2 = kx*kx + ky*ky + kz*kz;
          if (k_2 >= kmin*kmin and k_2 <= kmax*kmax) {
            nForcingModes ++;
            k2Disovec.push_back({kx, ky, kz});
            modeNames.push_back({"I" + std::to_string(nx) + std::to_string(ny) + std::to_string(nz), "J" + std::to_string(nx) + std::to_string(ny) + std::to_string(nz), "K" + std::to_string(nx) + std::to_string(ny) + std::to_string(nz)});
          }
        }
      }
    }
    k2DisoHost = IdefixHostArray2D<real>("k2DisoHost", nForcingModes, 3);
    k2Diso = IdefixArray2D<real>("k2Diso", nForcingModes, 3);
    for (int l=0; l<nForcingModes; l++) {
      k2DisoHost(l,0) = k2Disovec[l][0];
      k2DisoHost(l,1) = k2Disovec[l][1];
      k2DisoHost(l,2) = k2Disovec[l][2];
    }
    Kokkos::deep_copy(k2Diso, k2DisoHost);
    this->forcingModesIdir = IdefixArray4D<Kokkos::complex<real>>("forcingModesIdir", nForcingModes,
                              data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);
    this->forcingModesJdir = IdefixArray4D<Kokkos::complex<real>>("forcingModesJdir", nForcingModes,
                              data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);
    this->forcingModesKdir = IdefixArray4D<Kokkos::complex<real>>("forcingModesKdir", nForcingModes,
                              data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);
  }

  else if (input.CheckEntry("Forcing", "vsh")>=0) {
    this->forcingType = vsh;
    #if GEOMETRY != SPHERICAL
      IDEFIX_ERROR("You cannot have vsh forcing in another geometry than spherical.");
    #endif //GEOMETRY != SPHERICAL
    #if VSH == NO
      IDEFIX_ERROR("You cannot have vsh forcing without the VSH module.");
    #endif //VSH == NO
// if () IDEFIX_ERROR("You cannot have vsh forcing without the full sphere.");
    this->ellmin = input.Get<int>("Forcing","vsh", 0);
    this->ellmax = input.Get<int>("Forcing","vsh", 1);
    this->mmin = input.GetOrSet<int>("Forcing","vsh", 2, 0);
    this->mmax = input.GetOrSet<int>("Forcing","vsh", 3, this->ellmax);

    std::vector<std::vector<int>> ellmVshvec;
    for (int ell=ellmin; ell<ellmax; ell++) {
      for (int m=mmin; m<mmax & m<ell+1 ; m++) {
        this->nForcingModes++;
        ellmVshvec.push_back({ell,m});
        modeNames.push_back({"Y" + std::to_string(ell) + std::to_string(m), "S" + std::to_string(ell) + std::to_string(m), "T" + std::to_string(ell) + std::to_string(m)});
      }
    }
    ellmVshHost = IdefixHostArray2D<int>("ellmVshHost", nForcingModes, 2);
    ellmVsh = IdefixArray2D<int>("ellmVsh", nForcingModes, 2);
    for (int l=0; l<nForcingModes; l++) {
      ellmVshHost(l,0) = ellmVshvec[l][0];
      ellmVshHost(l,1) = ellmVshvec[l][1];
    }
    Kokkos::deep_copy(ellmVsh, ellmVshHost);

    this->forcingModesIdir = IdefixArray4D<Kokkos::complex<real>>("forcingModesIdir", nForcingModes,
                              data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);
    // WARNING, vsh case needs a special treatment since the forcing modes between
    // the theta and phi directions are coupled to each other...
    this->forcingModesJdir = IdefixArray4D<Kokkos::complex<real>>("forcingModesJdir", 2*nForcingModes,
                              data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);
    this->forcingModesKdir = IdefixArray4D<Kokkos::complex<real>>("forcingModesKdir", 2*nForcingModes,
                              data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);
  }

//  else { this->haveUserDef = true;}
  else { this->forcingType = userDef;} // TO BE CODED THOUGH

  if (input.CheckEntry("Forcing", "solenoidal")>=0) {
    #if GEOMETRY == POLAR or GEOMETRY == CYLINDRICAL
      IDEFIX_ERROR("Cannot have solenoidal forcing in polar in cylindrical geometry.");
    #endif //GEOMETRY == POLAR or GEOMETRY == CYLINDRICAL
    if (this->forcingType == vsh) IDEFIX_ERROR("Cannot have solenoidal forcing with vsh forcing. Select yourself the toroidal vector spherical harmonics in this case.");
    else if (this->forcingType == userDef) IDEFIX_ERROR("Cannot have solenoidal forcing with userdef forcing.");
    else this->haveSolenoidalForcing = true;
  }

  // Allocate required arrays
  this->forcingTerm = IdefixArray4D<real>("forcingTerm", COMPONENTS,
                              data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);
  this->pristineForcingTerm = IdefixArray4D<real>("pristineForcingTerm", COMPONENTS,
                              data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);
  this->solenoidalForcingTerm = IdefixArray4D<real>("solenoidalForcingTerm", COMPONENTS,
                              data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);
//  this->compressiveForcingTerm = IdefixArray4D<real>("compressiveForcingTerm", COMPONENTS,
//                              data->np_tot[KDIR], data->np_tot[JDIR], data->np_tot[IDIR]);
  this->means = IdefixArray2D<real>("means", nForcingModes, COMPONENTS);
  this->tcorrs = IdefixArray2D<real>("tcorrs", nForcingModes, COMPONENTS);
  this->epsilons = IdefixArray2D<real>("epsilons", nForcingModes, COMPONENTS);

// WARNING DO SOMETHING WITH MEANS TCORRS AND EPSILONS
  IdefixHostArray2D<real> mean("mean", nForcingModes, COMPONENTS);
  IdefixHostArray2D<real> tcorr("tcorr", nForcingModes, COMPONENTS);
  IdefixHostArray2D<real> epsilon("epsilon", nForcingModes, COMPONENTS);
  for (int l=0; l<nForcingModes; l++) {
    for (int dir=IDIR; dir<COMPONENTS; dir++) {
      mean(l,dir) = ZERO_F;
      tcorr(l,dir) = input.Get<real>("Forcing", "t_corr", 0);
      epsilon(l,dir) = input.Get<real>("Forcing", "epsilon", 0);
    }
  }
  Kokkos::deep_copy(this->means, mean);
  Kokkos::deep_copy(this->tcorrs, tcorr);
  Kokkos::deep_copy(this->epsilons, epsilon);
  this->oUprocesses.InitProcesses(folder, this->seed, this->nForcingModes, this->modeNames, this->means, this->tcorrs, this->epsilons);

  idfx::popRegion();
}

void Forcing::ShowConfig() {
  idfx::cout << "Forcing: ENABLED with seed " << seed << "." << std::endl;
  switch(forcingType) {
    case ForcingType::iso3D:
      idfx::cout << "Forcing: 3D isotropic." << std::endl;
      idfx::cout << "Forcing: kmin=" << kmin << " and kmax=" << kmax << " ." << std::endl;
      if (haveSolenoidalForcing) idfx::cout << "Forcing: solenoidal." << std::endl;
      break;
    case ForcingType::iso2D:
      idfx::cout << "Forcing: 2D isotropic with normal " << normal2DisoStr << "." << std::endl;
      idfx::cout << "Forcing: kmin=" << kmin << " and kmax=" << kmax << " ." << std::endl;
      if (haveSolenoidalForcing) idfx::cout << "Forcing: solenoidal." << std::endl;
      break;
    case ForcingType::vsh:
      idfx::cout << "Forcing: vector spherical harmonics." << std::endl;
      idfx::cout << "Forcing: ellmin=" << ellmin << " and ellmax=" << ellmax << " ." << std::endl;
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
//  IdefixArray4D<Kokkos::complex<real>> forcingModes = this->forcingModes;
  IdefixArray4D<Kokkos::complex<real>> forcingModesIdir = this->forcingModesIdir;
  IdefixArray4D<Kokkos::complex<real>> forcingModesJdir = this->forcingModesJdir;
  IdefixArray4D<Kokkos::complex<real>> forcingModesKdir = this->forcingModesKdir;
  IdefixArray1D<real> x1 = data->x[IDIR];
  IdefixArray1D<real> x2 = data->x[JDIR];
  IdefixArray1D<real> x3 = data->x[KDIR];
  IdefixArray2D<real> k3Diso = this->k3Diso;
  IdefixArray2D<real> k2Diso = this->k2Diso;
  IdefixArray2D<int> ellmVsh = this->ellmVsh;
  Kokkos::complex unit_j(0.,1.);
  int nForcingModes = this->nForcingModes;
  ForcingType forcingType = this->forcingType;
  idefix_for("Forcing::InitForcingModes",
              0, nForcingModes,
              0, data->np_tot[KDIR],
              0, data->np_tot[JDIR],
              0, data->np_tot[IDIR],
              KOKKOS_LAMBDA(int l, int k, int j, int i) {
                forcingModesIdir(l,k,j,i) = ZERO_F;
                forcingModesJdir(l,k,j,i) = ZERO_F;
                forcingModesKdir(l,k,j,i) = ZERO_F;
                if (forcingType == vsh) {
                  forcingModesJdir(l+nForcingModes,k,j,i) = ZERO_F;
                  forcingModesKdir(l+nForcingModes,k,j,i) = ZERO_F;
                }
              });
  switch(forcingType) {
    case ForcingType::iso3D:
      idefix_for("iso3D", 0, nForcingModes, 0, data->np_tot[KDIR], 0, data->np_tot[JDIR], 0, data->np_tot[IDIR],
                  KOKKOS_LAMBDA (int l, int k, int j, int i) {
                    real kx = k3Diso(l, IDIR);
                    real ky = k3Diso(l, JDIR);
                    real kz = k3Diso(l, KDIR);
                    forcingModesIdir(l,k,j,i) = exp(unit_j*(kx*x1(i) + ky*x2(j) + kz*x3(k)));
                    forcingModesJdir(l,k,j,i) = exp(unit_j*(kx*x1(i) + ky*x2(j) + kz*x3(k)));
                    forcingModesKdir(l,k,j,i) = exp(unit_j*(kx*x1(i) + ky*x2(j) + kz*x3(k)));
      });
//      Kokkos::deep_copy(forcingModesIdir, forcingModes);
//      Kokkos::deep_copy(forcingModesJdir, forcingModes);
//      Kokkos::deep_copy(forcingModesKdir, forcingModes);
      break;
    case ForcingType::iso2D:
      idefix_for("iso2D", 0, nForcingModes, 0, data->np_tot[KDIR], 0, data->np_tot[JDIR], 0, data->np_tot[IDIR],
                  KOKKOS_LAMBDA (int l, int k, int j, int i) {
                    real kx = k2Diso(l, IDIR);
                    real ky = k2Diso(l, JDIR);
                    real kz = k2Diso(l, KDIR);
                    forcingModesIdir(l,k,j,i) = exp(unit_j*(kx*x1(i) + ky*x2(j) + kz*x3(k)));
                    forcingModesJdir(l,k,j,i) = exp(unit_j*(kx*x1(i) + ky*x2(j) + kz*x3(k)));
                    forcingModesKdir(l,k,j,i) = exp(unit_j*(kx*x1(i) + ky*x2(j) + kz*x3(k)));
      });
//      Kokkos::deep_copy(forcingModesIdir, forcingModes);
//      Kokkos::deep_copy(forcingModesJdir, forcingModes);
//      Kokkos::deep_copy(forcingModesKdir, forcingModes);
      break;
    #if VSH == YES
      case ForcingType::vsh:
        IdefixArray4D<Kokkos::complex<real>> Ylm_r = this->data->Ylm_r;
        IdefixArray4D<Kokkos::complex<real>> Slm_th = this->data->Slm_th;
        IdefixArray4D<Kokkos::complex<real>> Slm_phi = this->data->Slm_phi;
        IdefixArray4D<Kokkos::complex<real>> Tlm_th = this->data->Tlm_th;
        IdefixArray4D<Kokkos::complex<real>> Tlm_phi = this->data->Tlm_phi;
        int nForcingModes = this->nForcingModes;
        idefix_for("vsh", 0, nForcingModes, 0, data->np_tot[KDIR], 0, data->np_tot[JDIR], 0, data->np_tot[IDIR],
                    KOKKOS_LAMBDA (int l, int k, int j, int i) {
                      int ell = ellmVsh(l, 0);
                      int m = ellmVsh(l, 1);
                      forcingModesIdir(l,k,j,i) = Ylm_r(ell,m,k,j);
                      forcingModesJdir(l,k,j,i) = Slm_th(ell,m,k,j);
                      forcingModesJdir(l+nForcingModes,k,j,i) = Tlm_th(ell,m,k,j);
                      forcingModesKdir(l,k,j,i) = Slm_phi(ell,m,k,j);
                      forcingModesKdir(l+nForcingModes,k,j,i) = Tlm_phi(ell,m,k,j);
        });
        break;
    #endif //VSH == YES
//    case ForcingType::userDef:
//      break;
  }
}

// This function compute the required forcing field
void Forcing::ComputeForcing(real dt) {
  idfx::pushRegion("Forcing::ComputeForcing");
  
  ResetForcingTerms();
  if (haveSolenoidalForcing) {
    ComputeSolenoidalForcing(dt);
    forcingTerm = solenoidalForcingTerm;
  } else {
    ComputePristineForcing(dt);
    forcingTerm = pristineForcingTerm;
  }
  oUprocesses.UpdateProcessesValues(dt);
  idfx::popRegion();
}

// This function compute the pristine forcing field
void Forcing::ComputePristineForcing(real dt) {
  idfx::pushRegion("Forcing::ComputePristineForcing");
  
  IdefixArray4D<Kokkos::complex<real>> forcingModesIdir = this->forcingModesIdir;
  IdefixArray4D<Kokkos::complex<real>> forcingModesJdir = this->forcingModesJdir;
  IdefixArray4D<Kokkos::complex<real>> forcingModesKdir = this->forcingModesKdir;
  IdefixArray4D<real> forcingTerm = this->forcingTerm;
  IdefixArray2D<Kokkos::complex<real>> ouValues = this->oUprocesses.ouValues;
  int nForcingModes = this->nForcingModes;
  idefix_for("ComputePristineForcing", 0, data->np_tot[KDIR], 0, data->np_tot[JDIR], 0, data->np_tot[IDIR],
              KOKKOS_LAMBDA (int k, int j, int i) {
                  Kokkos::complex<real> forcingX1=0.;
                  Kokkos::complex<real> forcingX2=0.;
                  Kokkos::complex<real> forcingX3=0.;
                  for (int l=0; l<nForcingModes; l++) {
                    forcingX1 += ouValues(l,IDIR)*forcingModesIdir(l,k,j,i);
                    #if COMPONENTS >= 2
                      forcingX2 += ouValues(l,JDIR)*forcingModesJdir(l,k,j,i);
                      #if VSH == YES
                        //in this case (l,JDIR) -> (l,Slm) and (l,KDIR) -> (l,Tlm)
                        forcingX2 += ouValues(l,KDIR)*forcingModesJdir(l+nForcingModes,k,j,i);
                      #endif //VSH == YES
                      #if COMPONENTS == 3
                        #if VSH == YES
                          //don't forget (l,JDIR) -> (l,Slm) and (l,KDIR) -> (l,Tlm)
                          forcingX3 += ouValues(l,JDIR)*forcingModesKdir(l,k,j,i);
                          forcingX3 += ouValues(l,KDIR)*forcingModesJdir(l+nForcingModes,k,j,i);
                        #else
                          forcingX3 += ouValues(l,KDIR)*forcingModesKdir(l,k,j,i);
                        #endif //VSH == YES
                      #endif //COMPONENTS == 3
                    #endif //COMPONENTS >= 2
                  }
                  forcingTerm(IDIR,k,j,i) += forcingX1.real();
                  forcingTerm(JDIR,k,j,i) += forcingX2.real();
                  forcingTerm(KDIR,k,j,i) += forcingX3.real();
  });

  idfx::popRegion();
}

// This function compute the solenoidal part of the forcing field
void Forcing::ComputeSolenoidalForcing(real dt) {
  idfx::pushRegion("Forcing::ComputeSolenoidalForcing");
  
  IdefixArray2D<real> kiso = (this->forcingType == iso3D) ? this->k3Diso : this->k2Diso;
  IdefixArray4D<Kokkos::complex<real>> forcingModesIdir = this->forcingModesIdir;
  IdefixArray4D<Kokkos::complex<real>> forcingModesJdir = this->forcingModesJdir;
  IdefixArray4D<Kokkos::complex<real>> forcingModesKdir = this->forcingModesKdir;
  IdefixArray4D<real> solenoidalForcingTerm = this->solenoidalForcingTerm;
  IdefixArray2D<Kokkos::complex<real>> ouValues = this->oUprocesses.ouValues;
  int nForcingModes = this->nForcingModes;
  Kokkos::complex unit_j(0.,1.);
  #if GEOMETRY == SPHERICAL
    IdefixArray1D<real> x1 = data->x[IDIR];
    IdefixArray1D<real> sinx2 = data->sinx2;
  #endif //GEOMETRY == SPHERICAL
  idefix_for("ComputeSolenoidalForcing", 0, data->np_tot[KDIR], 0, data->np_tot[JDIR], 0, data->np_tot[IDIR],
              KOKKOS_LAMBDA (int k, int j, int i) {
                  #if GEOMETRY == SPHERICAL
                    real r = x1(i);
                    real sin_1 = 1./sinx2(j);
                  #endif //GEOMETRY == SPHERICAL
                  Kokkos::complex<real> forcingX1 = 0.;
                  Kokkos::complex<real> forcingX2 = 0.;
                  Kokkos::complex<real> forcingX3 = 0.;
                  real kx1, kx2, kx3, sqrt_kk;
                  for (int l=0; l<nForcingModes; l++) {
                    kx1 = kiso(l,0);
                    kx2 = kiso(l,1);
                    kx3 = kiso(l,2);
                    #if GEOMETRY == CARTESIAN
                      sqrt_kk = sqrt(kx1*kx1 + kx2*kx2 + kx3*kx3);
                      //Here we are computing unit k
                      kx1 /= sqrt_kk;
                      kx2 /= sqrt_kk;
                      kx3 /= sqrt_kk;
                      //Here we are coding, for each mode k, a_k x ik/k
                      forcingX1 += unit_j*(ouValues(l,JDIR)*kx3 - ouValues(l,KDIR)*kx2)*forcingModesIdir(l,k,j,i);
                      #if COMPONENTS >= 2
                        forcingX2 += unit_j*(ouValues(l,KDIR)*kx1 - ouValues(l,IDIR)*kx3)*forcingModesJdir(l,k,j,i);
                        #if COMPONENTS == 3
                          forcingX3 += unit_j*(ouValues(l,IDIR)*kx2 - ouValues(l,JDIR)*kx1)*forcingModesKdir(l,k,j,i);
                        #endif //COMPONENTS == 3
                      #endif //COMPONENTS >= 2
                    #elif GEOMETRY == SPHERICAL
                      sqrt_kk = sqrt(kx1*kx1 + 1./(r*r)*(kx2*kx2 + sin_1*kx3*sin_1*kx3));
                      //Here we are computing unit k
                      kx1 /= sqrt_kk;
                      kx2 /= sqrt_kk;
                      kx3 /= sqrt_kk;
                      //Here we are coding, for each mode k, a_k x ik/k
                      forcingX1 += unit_j*(ouValues(l,JDIR)*sin_1/r*kx3 - ouValues(l,KDIR)/r*kx2)*forcingModesIdir(l,k,j,i);
                      #if COMPONENTS >= 2
                        forcingX2 += unit_j*(ouValues(l,KDIR)*kx1 - ouValues(l,IDIR)*sin_1/r*kx3)*forcingModesJdir(l,k,j,i);
                        #if COMPONENTS == 3
                          forcingX3 += unit_j*(ouValues(l,IDIR)/r*kx2 - ouValues(l,JDIR)*kx1)*forcingModesKdir(l,k,j,i);
                        #endif //COMPONENTS == 3
                      #endif //COMPONENTS >= 2
                    #endif //GEOMETRY == SPHERICAL
                    }
                  solenoidalForcingTerm(IDIR,k,j,i) += forcingX1.real();
                  solenoidalForcingTerm(JDIR,k,j,i) += forcingX2.real();
                  solenoidalForcingTerm(KDIR,k,j,i) += forcingX3.real();
  });

  idfx::popRegion();
}

//// This function compute the compressive part of the forcing field
//void Forcing::ComputeCompressiveForcing(real dt) {
//  idfx::pushRegion("Forcing::ComputeCompressiveForcing");
//  
//  ResetForcingTerms();
//
//  IdefixArray2D<real> kiso = (this->forcingType == iso3D) ? this->k3Diso : this->k2Diso;
//  IdefixArray4D<Kokkos::complex<real>> forcingModesIdir = this->forcingModesIdir;
//  IdefixArray4D<Kokkos::complex<real>> forcingModesJdir = this->forcingModesJdir;
//  IdefixArray4D<Kokkos::complex<real>> forcingModesKdir = this->forcingModesKdir;
//  IdefixArray4D<real> compressiveForcingTerm = this->compressiveForcingTerm;
//  IdefixArray2D<Kokkos::complex<real>> ouValues = this->oUprocesses.ouValues;
//  int nForcingModes = this->nForcingModes;
//  Kokkos::complex unit_j(0.,1.);
//  idefix_for("ComputeCompressiveForcing", 0, data->np_tot[KDIR], 0, data->np_tot[JDIR], 0, data->np_tot[IDIR],
//              KOKKOS_LAMBDA (int k, int j, int i) {
//                  Kokkos::complex<real> forcingX1 = 0.;
//                  Kokkos::complex<real> forcingX2 = 0.;
//                  Kokkos::complex<real> forcingX3 = 0.;
//                  real kx1, kx2, kx3;
//                  for (int l=0; l<nForcingModes; l++) {
//                    kx1 = kiso(l,0);
//                    kx2 = kiso(l,1);
//                    kx3 = kiso(l,2);
//                    kk = kx1*kx1 + kx2*kx2 + kx3*kx3;
//                    //Here we are computing unit k
//                    kx1 /= kk;
//                    kx2 /= kk;
//                    kx3 /= kk;
//                    //Here we are coding, for each mode k, a_k . ik/k
//                    forcingX1 += unit_j*kx1*ouValues(l,IDIR)*forcingModesIdir(l,k,j,i);
//                    #if COMPONENTS >= 2
//                      forcingX2 += unit_j*(ouValues(l,KDIR)*kx1 - ouValues(l,IDIR)*kx3)*forcingModesJdir(l,k,j,i);
//                      #if COMPONENTS == 3
//                        forcingX3 += unit_j*(ouValues(l,IDIR)*kx2 - ouValues(l,JDIR)*kx1)*forcingModesKdir(l,k,j,i);
//                      #endif //COMPONENTS == 3
//                    #endif //COMPONENTS >= 2
//                  }
//                  compressiveForcingTerm(IDIR,k,j,i) += forcingX1.real();
//                  compressiveForcingTerm(JDIR,k,j,i) += forcingX2.real();
//                  compressiveForcingTerm(KDIR,k,j,i) += forcingX3.real();
//  });
//
//  idfx::popRegion();
//}

// Fill the forcing term with zeros
void Forcing::ResetForcingTerms() {
  idfx::pushRegion("Forcing::ResetForcingTerms");
  IdefixArray4D<real> forcingTerm = this->forcingTerm;
  IdefixArray4D<real> solenoidalForcingTerm = this->solenoidalForcingTerm;
//  IdefixArray4D<real> compressiveForcingTerm = this->compressiveForcingTerm;
  idefix_for("Forcing::ResetForcingTerms",
              0, data->np_tot[KDIR],
              0, data->np_tot[JDIR],
              0, data->np_tot[IDIR],
              KOKKOS_LAMBDA(int k, int j, int i) {
                forcingTerm(IDIR,k,j,i) = ZERO_F;
                forcingTerm(JDIR,k,j,i) = ZERO_F;
                forcingTerm(KDIR,k,j,i) = ZERO_F;
                solenoidalForcingTerm(IDIR,k,j,i) = ZERO_F;
                solenoidalForcingTerm(JDIR,k,j,i) = ZERO_F;
                solenoidalForcingTerm(KDIR,k,j,i) = ZERO_F;
//                compressiveForcingTerm(IDIR,k,j,i) = ZERO_F;
//                compressiveForcingTerm(JDIR,k,j,i) = ZERO_F;
//                compressiveForcingTerm(KDIR,k,j,i) = ZERO_F;
              });
  idfx::popRegion();
}

