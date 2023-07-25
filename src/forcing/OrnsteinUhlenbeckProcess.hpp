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

#include <iostream>
#include <random>
#include "idefix.hpp"

class OrnsteinUhlenbeckProcess {
private:
    real mean;
    real tcorr;
//    real sigma;
//    real dt;
    real ou;

    std::default_random_engine generator;
    std::normal_distribution<real> normal_distribution;

public:
    OrnsteinUhlenbeckProcess(); // Default (empty) constructor
//    InitProcess(int, real, real, real, real);
    void InitProcess(int, real, real);
//    real getNextValue();
    real GetNextValue(real, real);
};

class OrnsteinUhlenbeckProcesses {
private:
    int lmax;
    int mmax;

    real mean;
    real tcorr;
//    real sigma;
//    real dt;

    std::vector<OrnsteinUhlenbeckProcess> ouProcesses;
public:
    IdefixArray3D<real> ouValues;

    OrnsteinUhlenbeckProcesses(); // Default (empty) constructor
//    InitProcess(int, real, real, real, real);
    void InitProcesses(int, int, real, real);
//    real getNextValue();
    void UpdateProcesses(real, real);
};

#endif  // ORNSTEIN_UHLENBECK_PROCESS_HPP
