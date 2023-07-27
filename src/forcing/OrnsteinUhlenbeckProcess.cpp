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

OrnsteinUhlenbeckProcess::OrnsteinUhlenbeckProcess()
{ // Default (empty) constructor
}

//void OrnsteinUhlenbeckProcess::InitProcess(int seed, real mean, real tcorr, real sigma, real dt) {
void OrnsteinUhlenbeckProcess::InitProcess(int seed, real mean, real tcorr) {
  this->mean = mean;
  this->tcorr = tcorr;
  this->ou = mean;
  this->generator = std::default_random_engine(seed);
  this->normal_distribution = std::normal_distribution(0.,1.);
}

//real OrnsteinUhlenbeckProcess::GetNextValue() {
real OrnsteinUhlenbeckProcess::GetNextValue(real dt, real sigma) {
  real normal = this->normal_distribution(this->generator);
  real expTerm = std::exp(-dt/this->tcorr);
  real dou = sigma*std::sqrt((1 - expTerm*expTerm)*this->tcorr/2.)*normal;
  ou = ou*expTerm + dou;
  return ou;
}

OrnsteinUhlenbeckProcesses::OrnsteinUhlenbeckProcesses()
{ // Default (empty) constructor
}

//void OrnsteinUhlenbeckProcess::InitProcess(int seed, real mean, real tcorr, real sigma, real dt) {
void OrnsteinUhlenbeckProcesses::InitProcesses(int lmax, int mmax, real mean, real tcorr) {
  this->lmax = lmax;
  this->mmax = mmax;
  this->mean = mean;
  this->tcorr = tcorr;
  this->ouProcesses = std::vector<OrnsteinUhlenbeckProcess>(3*lmax*mmax);
  this->ouValues = IdefixArray3D<real> ("ouProcesses", 3, lmax, mmax);

  for (int comp=0; comp<3; comp++) {
    for (int l=0; l<lmax; l++) {
      for (int m=0; m<mmax; m++) {
        int seed = pow(2,l)*(2*m+1);
        ouProcesses[comp*(lmax*mmax)+l*mmax+m].InitProcess(seed,mean,tcorr);
        ouValues(comp,l,m) = mean;
      }
    }
  }
}

void OrnsteinUhlenbeckProcesses::UpdateProcesses(real dt, real epsilon_Ylm, real epsilon_Slm, real epsilon_Tlm) {
  real value;
  for (int l=0; l<this->lmax; l++) {
    for (int m=0; m<this->mmax; m++) {
      value = ouProcesses[0*(this->lmax*this->mmax)+l*this->mmax+m].GetNextValue(dt, epsilon_Ylm);
      ouValues(0,l,m) = value;
      value = ouProcesses[1*(this->lmax*this->mmax)+l*this->mmax+m].GetNextValue(dt, epsilon_Slm);
      ouValues(1,l,m) = value;
      value = ouProcesses[2*(this->lmax*this->mmax)+l*this->mmax+m].GetNextValue(dt, epsilon_Tlm);
      ouValues(2,l,m) = value;
    }
  }
//        std::cout << "Done" << std::endl;
}

