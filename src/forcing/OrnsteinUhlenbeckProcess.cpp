// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2022 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2021 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#include "OrnsteinUhlenbeckProcess.hpp"

OrnsteinUhlenbeckProcesses::OrnsteinUhlenbeckProcesses()
{ // Default (empty) constructor
}

void OrnsteinUhlenbeckProcesses::InitProcesses(int seed, int lmin, int lmax, int mmin, int mmax, real mean, real tcorr, real eps_Ylm, real eps_Slm, real eps_Tlm) {
  this->lmin = lmin;
  this->lmax = lmax;
  this->mmin = mmin;
  this->mmax = mmax;
  this->epsilons = IdefixArray3D<real> ("ouEpsilons", 3,lmax,mmax);
  this->tcorrs = IdefixArray3D<real> ("ouTcorrs", 3,lmax,mmax);
  this->means = IdefixArray3D<real> ("ouMeans", 3,lmax,mmax);
  this->ouValues = IdefixArray3D<real> ("ouProcesses", 3, lmax, mmax);
  this->random_pool = Kokkos::Random_XorShift64_Pool<> (/*seed=*/seed);

  IdefixArray3D<real> means = this->means;
  IdefixArray3D<real> tcorrs = this->tcorrs;
  IdefixArray3D<real> epsilons = this->epsilons;
  IdefixArray3D<real> ouValues = this->ouValues;
  idefix_for("InitProcesses", 0, 3, 0, lmax, 0, mmax,
              KOKKOS_LAMBDA (int vshcomp, int l, int m) {
        if (l >= lmin & m >= mmin) {
          means(vshcomp,l,m) = mean;
          tcorrs(vshcomp,l,m) = tcorr;
          if (vshcomp==0) { epsilons(vshcomp,l,m) = eps_Ylm;}
          if (vshcomp==1) { epsilons(vshcomp,l,m) = eps_Slm;}
          if (vshcomp==2) { epsilons(vshcomp,l,m) = eps_Tlm;}
          ouValues(vshcomp,l,m) = mean;
        } else {
          means(vshcomp,l,m) = ZERO_F;
          tcorrs(vshcomp,l,m) = ZERO_F;
          epsilons(vshcomp,l,m) = ZERO_F;
        }
  });
}

//void OrnsteinUhlenbeckProcesses::UpdateProcessesValues(real dt, IdefixArray1D<real> epsilons) {
void OrnsteinUhlenbeckProcesses::UpdateProcessesValues(real dt) {
  IdefixArray3D<real> means = this->means;
  IdefixArray3D<real> tcorrs = this->tcorrs;
  IdefixArray3D<real> epsilons = this->epsilons;
  IdefixArray3D<real> ouValues = this->ouValues;
  Kokkos::Random_XorShift64_Pool<> random_pool = this->random_pool;
  idefix_for("UpdateProcesses", 0, 3, lmin, lmax, mmin, mmax,
              KOKKOS_LAMBDA (int vshcomp, int l, int m) {
      if (m < l+1) {
        auto generator = random_pool.get_state();
        real normal = generator.normal(0., 1.);
  printf("%f;", normal);
        random_pool.free_state(generator);
        real expTerm = std::exp(-dt/tcorrs(vshcomp,l,m));
        real dou = std::sqrt(epsilons(vshcomp,l,m)/tcorrs(vshcomp,l,m)*(1 - expTerm*expTerm))*normal;
        real newValue = means(vshcomp,l,m) + (ouValues(vshcomp,l,m)-means(vshcomp,l,m))*expTerm + dou;
        ouValues(vshcomp,l,m) = newValue;
  std::string ouFile = "checkFiles/ou_comp" + std::to_string(vshcomp) + "_lmax" + std::to_string(l) + "_mmax" + std::to_string(m) + ".csv";
  myfile.open (ouFile, std::fstream::app);
  myfile << newValue << ";";
  myfile.close();
  std::string normalFile = "checkFiles/normal_comp" + std::to_string(vshcomp) + "_lmax" + std::to_string(l) + "_mmax" + std::to_string(m) + ".csv";
  myfile.open (normalFile, std::fstream::app);
  myfile << normal << ";";
  myfile.close();
      } else {
        ouValues(vshcomp,l,m) = 0.;
      }
  });
}

