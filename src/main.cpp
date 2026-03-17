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

#include "myIdefix.hpp"

int main( int argc, char* argv[] ) {
  MyIdefix idefix(argc, argv);
  idefix.Initialise(argc, argv);
  idefix.DoMainLoop();
  idefix.Finalise();
  return(0);
}
