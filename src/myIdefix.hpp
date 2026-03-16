// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef MYIDEFIX_HPP_
#define MYIDEFIX_HPP_
#include "idefix.hpp"
#include "input.hpp"
#include "grid.hpp"
#include "gridHost.hpp"
#include "dataBlock.hpp"
#include "timeIntegrator.hpp"
#include "setup.hpp"
#ifdef WITH_PYTHON
#include "pydefix.hpp"
#endif
#include "output.hpp"

class MyIdefix {
 public:
  MyIdefix(int, char**);
  void Initialise(int, char**);
  void DoMainLoop();
  void Finalise();

 private:
  bool initKokkosBeforeMPI;
  int returnCode;
  real tstop;
  Kokkos::Timer timer;
  Input input;
  Grid grid;
  GridHost gridHost;
  DataBlock data;
  TimeIntegrator Tint;
#ifdef WITH_PYTHON
  Pydefix pydefix;
#endif
  Output output;
  Setup mysetup;
};
#endif // OUTPUT_OUTPUT_HPP_

