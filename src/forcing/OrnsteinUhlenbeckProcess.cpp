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

void OrnsteinUhlenbeckProcesses::InitProcesses(int lmin, int lmax, int mmin, int mmax, real mean, real tcorr, real eps_Ylm, real eps_Slm, real eps_Tlm) {
  this->lmin = lmin;
  this->lmax = lmax;
  this->mmin = mmin;
  this->mmax = mmax;
  this->epsilons = IdefixArray3D<real> ("ouEpsilons", 3,lmax,mmax);
  this->tcorrs = IdefixArray3D<real> ("ouTcorrs", 3,lmax,mmax);
  this->means = IdefixArray3D<real> ("ouMeans", 3,lmax,mmax);
  this->ouValues = IdefixArray3D<real> ("ouProcesses", 3, lmax, mmax);
  this->random_pool = Kokkos::Random_XorShift64_Pool<> (/*seed=*/1);

  IdefixArray3D<real> means = this->means;
  IdefixArray3D<real> tcorrs = this->tcorrs;
  IdefixArray3D<real> epsilons = this->epsilons;
  IdefixArray3D<real> ouValues = this->ouValues;
  idefix_for("InitProcesses", 0, 3, 0, lmax, 0, mmax,
              KOKKOS_LAMBDA (int vshcomp, int l, int m) {
//        int seed = pow(2,l)*(2*m+1);
//        this->generator = std::default_random_engine(seed);
//        this->normal_distribution = std::normal_distribution(0.,1.);
        means(vshcomp,l,m) = mean;
        tcorrs(vshcomp,l,m) = tcorr;
        if (vshcomp==0) { epsilons(vshcomp,l,m) = eps_Ylm;}
        if (vshcomp==1) { epsilons(vshcomp,l,m) = eps_Slm;}
        if (vshcomp==2) { epsilons(vshcomp,l,m) = eps_Tlm;}
        ouValues(vshcomp,l,m) = mean;
  });
}

//void OrnsteinUhlenbeckProcesses::UpdateProcessesValues(real dt, IdefixArray1D<real> epsilons) {
void OrnsteinUhlenbeckProcesses::UpdateProcessesValues(real dt) {
  IdefixArray3D<real> means = this->means;
  IdefixArray3D<real> tcorrs = this->tcorrs;
  IdefixArray3D<real> epsilons = this->epsilons;
  IdefixArray3D<real> ouValues = this->ouValues;
  Kokkos::Random_XorShift64_Pool<> random_pool = this->random_pool;
  idefix_for("UpdateProcesses", 0, 3, 0, lmax, 0, mmax,
              KOKKOS_LAMBDA (int vshcomp, int l, int m) {
//      real normal = 0.5;
   auto generator = random_pool.get_state();
   real normal = generator.normal(0., 1.);
   random_pool.free_state(generator);
      real expTerm = std::exp(-dt/tcorrs(vshcomp,l,m));
      real dou = std::sqrt(epsilons(vshcomp,l,m)/tcorrs(vshcomp,l,m)*(1 - expTerm*expTerm))*normal;
      real newValue = means(vshcomp,l,m) + (ouValues(vshcomp,l,m)-means(vshcomp,l,m))*expTerm + dou;
      ouValues(vshcomp,l,m) = newValue;
//      ouValues(vshcomp,l,m) = 10.;
  });
}

