// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

/*
//@HEADER
// ************************************************************************
//
//                        IDEFIX v 2.2.01
//
// ************************************************************************
//@HEADER
*/

#include <Kokkos_Core.hpp>

#include "menhir.hpp"

int main( int argc, char* argv[] ) {

  int initKokkosBeforeMPI = false;

  // return code is zero if the simulation reached final time
  // >0 if a fatal error occured (too small timestep, Nans)
  // <0 if simulation was interrupted (max_runtime or user-triggered interruption
  int returnCode = 0;

  // When running on GPUS with Omnipath network,
  // Kokkos needs to be initialised *before* the MPI layer
#ifdef KOKKOS_ENABLE_CUDA
  if(std::getenv("PSM2_CUDA") != NULL) {
    initKokkosBeforeMPI = true;
  }
#endif

  if(initKokkosBeforeMPI)  Kokkos::initialize( argc, argv );

#ifdef WITH_MPI
  MPI_Init(&argc,&argv);
#endif

  if(!initKokkosBeforeMPI) Kokkos::initialize( argc, argv );

{
  Menhir menhir(argc, argv);
  if (menhir.haveDNS) menhir.PerformDNS();
  else if (menhir.haveNewton) menhir.PerformNewton();
  else if (menhir.haveStability) menhir.PerformStability();
  else if (menhir.haveContinuation) menhir.PerformContinuation();
}

  if(returnCode<0) {
    idfx::cout << "Main: Job was interrupted before completion." << std::endl;
  } else if (returnCode>0) {
    idfx::cout << "Main: Job was aborted because of an unrecoverable error." << std::endl;
  } else {
    idfx::cout << "Main: Job completed successfully." << std::endl;
  }
  Kokkos::finalize();

#ifdef WITH_MPI
  MPI_Finalize();
#endif
  return(0);
}
