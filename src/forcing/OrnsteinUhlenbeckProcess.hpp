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
    int lmin;
    int mmin;
    int lmax;
    int mmax;

    IdefixArray3D<real> means;
    IdefixArray3D<real> tcorrs;
    IdefixArray3D<real> epsilons;

//    std::default_random_engine generator;
//    std::normal_distribution<real> normal_distribution;
    Kokkos::Random_XorShift64_Pool<> random_pool;

public:
    IdefixArray3D<real> ouValues;

    OrnsteinUhlenbeckProcesses(); // Default (empty) constructor
//    void InitProcesses(int, int, int, int, int, real, real, real, real, real);
    void InitProcesses(std::string, int, int, int, int, int, real, real, real, real, real);
    void UpdateProcessesValues(real);

    std::string filename;
    void ResetProcessesValues();
    void WriteProcessesValues();
    int precision;

    std::ofstream file;
};

#endif  // ORNSTEIN_UHLENBECK_PROCESS_HPP
