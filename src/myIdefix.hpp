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
  std::unique_ptr<Input> input;
  std::unique_ptr<Grid> grid;
  std::unique_ptr<GridHost> gridHost;
  std::unique_ptr<DataBlock> data;
  std::unique_ptr<TimeIntegrator> Tint;
#ifdef WITH_PYTHON
  std::unique_ptr<Pydefix> pydefix;
#endif
  std::unique_ptr<Output> output;
  std::unique_ptr<Setup> mysetup;
};
#endif // MYIDEFIX_HPP_

