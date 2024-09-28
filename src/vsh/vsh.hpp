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
    Vsh(DataBlockHost*, int); // Default constructor

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

    IdefixArray1D<real>::HostMirror x1; // spherical bessel function
    IdefixArray1D<real>::HostMirror x1l; // spherical bessel function

    IdefixArray2D<real>::HostMirror jl; // spherical bessel function
    IdefixArray2D<real>::HostMirror jls; // spherical bessel function at interface
    // WARNING: with shtns, Tlm equals -Tlm definition from wiki. See SHTNS website
    IdefixArray4D<Kokkos::complex<real>>::HostMirror Ylm_r;
    IdefixArray4D<Kokkos::complex<real>>::HostMirror Slm_th;
    IdefixArray4D<Kokkos::complex<real>>::HostMirror Slm_phi;
    IdefixArray4D<Kokkos::complex<real>>::HostMirror Tlm_th;
    IdefixArray4D<Kokkos::complex<real>>::HostMirror Tlm_phi;
    IdefixArray4D<Kokkos::complex<real>>::HostMirror Slm_ths;
    IdefixArray4D<Kokkos::complex<real>>::HostMirror Slm_phis;
    IdefixArray4D<Kokkos::complex<real>>::HostMirror Tlm_ths;
    IdefixArray4D<Kokkos::complex<real>>::HostMirror Tlm_phis;

    void ShowConfig();
    void Generatejl();
    void GenerateCellVsh();
    void GenerateInterfaceVsh();
    void GenerateCellGhostVsh();
    void GenerateInterfaceGhostVsh();
 
  private:
    int write;

    void write_vect(std::string, std::complex<real>*, int);
    void write_mx(std::string, std::complex<real>*, int, int);
};
#endif // VSH_VSH_HPP_
