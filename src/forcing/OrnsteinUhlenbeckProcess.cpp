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

void OrnsteinUhlenbeckProcesses::InitProcesses(std::string folder, int seed, int nForcingModes, IdefixArray1D<real> mean, IdefixArray1D<real> tcorr, IdefixArray1D<real> epsilon) {
  this->nForcingModes = nForcingModes;
  this->epsilons = IdefixArray1D<real> ("ouEpsilons", nForcingModes);
  this->tcorrs = IdefixArray1D<real> ("ouTcorrs", nForcingModes);
  this->means = IdefixArray1D<real> ("ouMeans", nForcingModes);
  this->ouValues = IdefixArray1D<real> ("ouValues", nForcingModes);
  this->normalValues = IdefixArray1D<real> ("normalValues", nForcingModes);
  this->random_pool = Kokkos::Random_XorShift64_Pool<> (/*seed=*/seed);

  IdefixArray1D<real> means = this->means;
  IdefixArray1D<real> tcorrs = this->tcorrs;
  IdefixArray1D<real> epsilons = this->epsilons;
  IdefixArray1D<real> ouValues = this->ouValues;
  IdefixArray1D<real> mean = mean;
  IdefixArray1D<real> tcorr = tcorr;
  IdefixArray1D<real> epsilon = epsilon;
  idefix_for("InitProcesses", 0, nForcingModes
              KOKKOS_LAMBDA (int l) {
        means(l) = mean(l);
        tcorrs(l) = tcorr(l);
        ouValues(l) = mean(l);
  });

  this->ouFilename = folder + "/ou_prank" + std::to_string(idfx::prank) + "_seed" + std::to_string(seed) + ".dat";
  this->normalFilename = folder + "/normal_prank" + std::to_string(idfx::prank) + "_seed" + std::to_string(seed) + ".dat";
  this->precision = 10;
  this->ouValuesHost = IdefixHostArray1D<real> ("ouValuesHost", l);
  this->normalValuesHost = IdefixHostArray1D<real> ("normalValuesHost", l);
}

//void OrnsteinUhlenbeckProcesses::UpdateProcessesValues(real dt, IdefixArray1D<real> epsilons) {
void OrnsteinUhlenbeckProcesses::UpdateProcessesValues(real dt) {
  IdefixArray1D<real> means = this->means;
  IdefixArray1D<real> tcorrs = this->tcorrs;
  IdefixArray1D<real> epsilons = this->epsilons;
  IdefixArray1D<real> ouValues = this->ouValues;
  IdefixArray1D<real> normalValues = this->normalValues;
  Kokkos::Random_XorShift64_Pool<> random_pool = this->random_pool;
  idefix_for("UpdateProcesses", 0, nForcingModes,
              KOKKOS_LAMBDA (int l) {
      auto generator = random_pool.get_state();
      real normal = generator.normal(0., 1.);
      random_pool.free_state(generator);
      normalValues(l) = normal;
      real expTerm = std::exp(-dt/tcorrs(l));
      real dou = std::sqrt(epsilons(l)/tcorrs(l)*(1. - expTerm*expTerm))*normal;
      real newValue = means(l) + (ouValues(l)-means(l))*expTerm + dou;
      ouValues(l) = newValue;
  });
}

void OrnsteinUhlenbeckProcesses::ResetProcessesValues() {
  if(idfx::prank==0) {
    file.open(ouFilename, std::ios::trunc);
    int col_width = precision + 10;
    file << std::setw(col_width) << "t";
    for (int l=0; l<nForcingModes; l++) {
      std::string current_name = std::to_string(l);
      file << std::setw(col_width) << current_name;
    }
    file << std::endl;
    file.close();
  }
}

void OrnsteinUhlenbeckProcesses::WriteProcessesValues(real time) {
  if(idfx::prank==0) {
    int col_width = precision + 10;
    file.open(ouFilename, std::ios::app);
    file.precision(precision);
    this->file << std::scientific << std::setw(col_width) << time;
    Kokkos::deep_copy(ouValuesHost, ouValues);
    for (int l=0; l<nForcingModes; l++) {
      file << std::scientific << std::setw(col_width) << ouValuesHost(l);
    }
    file << std::endl;
    file.close();
  }
}

void OrnsteinUhlenbeckProcesses::ResetNormalValues() {
  if(idfx::prank==0) {
    file.open(normalFilename, std::ios::trunc);
    int col_width = precision + 10;
    file << std::setw(col_width) << "t";
    for (int l=0; l<nForcingModes; l++) {
      std::string current_name = std::to_string(l);
      file << std::setw(col_width) << current_name;
    }
    file << std::endl;
    file.close();
  }
}

void OrnsteinUhlenbeckProcesses::WriteNormalValues(real time) {
  if(idfx::prank==0) {
    int col_width = precision + 10;
    file.open(normalFilename, std::ios::app);
    file.precision(precision);
    this->file << std::scientific << std::setw(col_width) << time;
    Kokkos::deep_copy(normalValuesHost, normalValues);
    for (int l=0; l<nForcingModes; l++) {
      file << std::scientific << std::setw(col_width) << normalValuesHost(l);
    }
    file << std::endl;
    file.close();
  }
}

