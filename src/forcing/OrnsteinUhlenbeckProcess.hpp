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
    real ou;

    std::default_random_engine generator;
    std::normal_distribution<real> normal_distribution;

public:
    OrnsteinUhlenbeckProcess(); // Default (empty) constructor
    void InitProcess(int, real, real);
    real GetNextValue(real, real);
};

class OrnsteinUhlenbeckProcesses {
private:
    int lmax;
    int mmax;

    real mean;
    real tcorr;

    std::vector<OrnsteinUhlenbeckProcess> ouProcesses;
public:
    IdefixArray3D<real> ouValues;

    OrnsteinUhlenbeckProcesses(); // Default (empty) constructor
    void InitProcesses(int, int, real, real);
    void UpdateProcesses(real, real, real, real);
};

#endif  // ORNSTEIN_UHLENBECK_PROCESS_HPP
