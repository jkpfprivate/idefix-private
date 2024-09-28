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

void OrnsteinUhlenbeckProcesses::InitProcesses(std::string folder, int seed, int nSeries, IdefixArray2D<real> mean, IdefixArray2D<real> tcorr, IdefixArray2D<real> epsilon) {
  this->nSeries = nSeries;
  this->epsilons = IdefixArray2D<real> ("ouEpsilons", nSeries, KDIR);
  this->tcorrs = IdefixArray2D<real> ("ouTcorrs", nSeries, KDIR);
  this->means = IdefixArray2D<real> ("ouMeans", nSeries, KDIR);
  this->ouValues = IdefixArray2D<Kokkos::complex<real>> ("ouValues", nSeries, KDIR);
  this->normalValues = IdefixArray2D<real> ("normalValues", nSeries, KDIR);
  this->random_pool = Kokkos::Random_XorShift64_Pool<> (/*seed=*/seed);

  IdefixArray2D<real> means = this->means;
  IdefixArray2D<real> tcorrs = this->tcorrs;
  IdefixArray2D<real> epsilons = this->epsilons;
  IdefixArray2D<Kokkos::complex<real>> ouValues = this->ouValues;
  idefix_for("InitProcesses", 0, nSeries, IDIR, KDIR,
              KOKKOS_LAMBDA (int l, int dir) {
        means(l, dir) = mean(l, dir);
        tcorrs(l, dir) = tcorr(l, dir);
        ouValues(l, dir) = mean(l, dir);
  });

  this->ouFilename = folder + "/ou_prank" + std::to_string(idfx::prank) + "_seed" + std::to_string(seed) + ".dat";
  this->normalFilename = folder + "/normal_prank" + std::to_string(idfx::prank) + "_seed" + std::to_string(seed) + ".dat";
  this->precision = 10;
  this->ouValuesHost = IdefixHostArray2D<Kokkos::complex<real>> ("ouValuesHost", nSeries);
  this->normalValuesHost = IdefixHostArray2D<real> ("normalValuesHost", nSeries);
}

//void OrnsteinUhlenbeckProcesses::UpdateProcessesValues(real dt, IdefixArray1D<real> epsilons) {
void OrnsteinUhlenbeckProcesses::UpdateProcessesValues(real dt) {
  IdefixArray2D<real> means = this->means;
  IdefixArray2D<real> tcorrs = this->tcorrs;
  IdefixArray2D<real> epsilons = this->epsilons;
  IdefixArray2D<Kokkos::complex<real>> ouValues = this->ouValues;
  IdefixArray2D<real> normalValues = this->normalValues;
  Kokkos::Random_XorShift64_Pool<> random_pool = this->random_pool;
  idefix_for("UpdateProcesses", 0, nSeries, IDIR, KDIR,
              KOKKOS_LAMBDA (int l, int dir) {
      auto generator = random_pool.get_state();
      real normal = generator.normal(0., 1.);
      random_pool.free_state(generator);
      normalValues(l, dir) = normal;
      real expTerm = std::exp(-dt/tcorrs(l, dir));
      real dou = std::sqrt(epsilons(l, dir)/tcorrs(l, dir)*(1. - expTerm*expTerm))*normal;
      real mod = std::sqrt(pow(ouValues(l,dir).real(), 2) + pow(ouValues(l,dir).imag(), 2));
      real newValue = means(l, dir) + (mod-means(l, dir))*expTerm + dou;

      real arg = (ouValues(l,dir).imag()==0. and ouValues(l,dir).real()<0.) ? M_PI : 2*atan(ouValues(l,dir).imag()/(ouValues(l,dir).real() + mod));
      dou = std::sqrt(2.*M_PI*(1. - expTerm*expTerm))*normal;
      real newArg = arg*expTerm + dou;
      Kokkos::complex<real> compNewValue(newValue*cos(newArg), newValue*sin(newArg));
      ouValues(l, dir) = compNewValue;
  });
}

void OrnsteinUhlenbeckProcesses::ResetProcessesValues() {
  if(idfx::prank==0) {
    file.open(ouFilename, std::ios::trunc);
    int col_width = precision + 10;
    file << std::setw(col_width) << "t";
    for (int l=0; l<nSeries; l++) {
      for (int dir=IDIR; dir<KDIR; dir++) {
        std::string current_name = std::to_string(l) + std::to_string(dir);
        file << std::setw(col_width) << current_name;
      }
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
    for (int l=0; l<nSeries; l++) {
      for (int dir=IDIR; dir<KDIR; dir++) {
        real mod = std::sqrt(pow(ouValues(l,dir).real(), 2) + pow(ouValues(l,dir).imag(), 2));
        file << std::scientific << std::setw(col_width) << mod;
      }
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
    for (int l=0; l<nSeries; l++) {
      for (int dir=IDIR; dir<KDIR; dir++) {
        std::string current_name = std::to_string(l) + std::to_string(dir);
        file << std::setw(col_width) << current_name;
      }
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
    for (int l=0; l<nSeries; l++) {
      for (int dir=IDIR; dir<KDIR; dir++) {
        file << std::scientific << std::setw(col_width) << normalValuesHost(l, dir);
      }
    }
    file << std::endl;
    file.close();
  }
}

