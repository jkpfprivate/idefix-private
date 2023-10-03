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

#ifndef VSH_VSH_HPP_
#define VSH_VSH_HPP_

#include "idefix.hpp"
//#include "input.hpp"
#include "shtns.h"

class DataBlockHost;

class Vsh {
  public:
//    Vsh(Input&, DataBlockHost*); // Default constructor
//    Vsh(int, int, DataBlockHost*); 
    Vsh(DataBlockHost*); 

    int nphi;
    int nphi_proc;
    int nphi_shtns;
    int ntheta;
    int ntheta_proc;
    int ntheta_shtns;
    int nr;
    int nr_proc;
    int koffset;
    int joffset;
    int kghost;
    int jghost;
    int ighost;
    int lmax;
    int mmax;
//    IdefixArray1D<real> x1; // spherical bessel function
//    IdefixArray1D<real> x1l; // spherical bessel function
    IdefixArray1D<real>::HostMirror x1; // spherical bessel function
    IdefixArray1D<real>::HostMirror x1l; // spherical bessel function

    IdefixArray2D<real> jl; // spherical bessel function
    IdefixArray2D<real> jls; // spherical bessel function at interface
    // WARNING: with shtns, Tlm equals -Tlm definition from wiki. See SHTNS website
    IdefixArray4D<real> Ylm_r, Slm_th, Slm_phi, Tlm_th, Tlm_phi;
    IdefixArray4D<real> Slm_ths, Slm_phis, Tlm_ths, Tlm_phis;

    void ShowConfig();
    void Generatejl();
    void GenerateCellVsh(int);
    void GenerateInterfaceVsh(int);
    void GenerateCellGhostVsh();
    void GenerateInterfaceGhostVsh();
 
  private:
    void write_vect(std::string, double *, int);
    void write_mx(std::string, double *, int, int);
};

#endif // VSH_VSH_HPP_
