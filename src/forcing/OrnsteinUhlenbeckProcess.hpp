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

#ifndef ORNSTEIN_UHLENBECK_PROCESS_HPP
#define ORNSTEIN_UHLENBECK_PROCESS_HPP

//#include <random>
//#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>
#include "idefix.hpp"
#include <iostream>
#include <fstream>
#include <iomanip>
//#include <ofstream>


class OrnsteinUhlenbeckProcesses {
private:
    int nForcingModes;

    IdefixArray1D<real> means;
    IdefixArray1D<real> tcorrs;
    IdefixArray1D<real> epsilons;

    Kokkos::Random_XorShift64_Pool<> random_pool;

public:
    IdefixArray1D<real> ouValues;
    IdefixHostArray1D<real> ouValuesHost;
    IdefixArray1D<real> normalValues;
    IdefixHostArray1D<real> normalValuesHost;

    OrnsteinUhlenbeckProcesses(); // Default (empty) constructor
//    void InitProcesses(int, int, int, int, int, real, real, real, real, real);
    void InitProcesses(std::string, int, int, IdefixArray1D<real>, IdefixArray1D<real>, IdefixArray1D<real>);
    void UpdateProcessesValues(real);

    std::string ouFilename;
    std::string normalFilename;
    void ResetProcessesValues();
    void WriteProcessesValues(real);
    void ResetNormalValues();
    void WriteNormalValues(real);
    int precision;

    std::ofstream file;
};

#endif  // ORNSTEIN_UHLENBECK_PROCESS_HPP
