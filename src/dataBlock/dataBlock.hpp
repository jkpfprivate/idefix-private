// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2022 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef DATABLOCK_DATABLOCK_HPP_
#define DATABLOCK_DATABLOCK_HPP_

#include <vector>
#include <string>
#include "idefix.hpp"
#include "grid.hpp"
#include "gridHost.hpp"
#include "hydro.hpp"
#include "fargo.hpp"
#include "gravity.hpp"
#include "stateContainer.hpp"

//TODO(lesurg) What is this standing for?
#define BOUNDARY_

class DataBlock {
 public:
  // Local grid information
  std::vector<IdefixArray1D<real>> x;    ///< geometrical central points
  std::vector<IdefixArray1D<real>> xr;   ///< cell right interface
  std::vector<IdefixArray1D<real>> xl;   ///< cell left interface
  std::vector<IdefixArray1D<real>> dx;   ///< cell width
  std::vector<IdefixArray1D<real>> xgc;  ///< cell geometrical cell center
  IdefixArray1D<real> rt;      ///< In spherical coordinates, gives $\tilde{r}$
  IdefixArray1D<real> sinx2m;  ///< In spherical coordinates,
                               ///< gives sin(th) at a j-1/2 interface
  IdefixArray1D<real> tanx2m;  ///< In spherical coordinates,
                               ///< gives tan(th) at a j-1/2 interface
  IdefixArray1D<real> sinx2;   ///< In spherical coordinates, gives sin(th) at the cell center
  IdefixArray1D<real> tanx2;   ///< In spherical coordinates, gives tan(th) at the cell center
  IdefixArray1D<real> dmu;     ///< In spherical coordinates,
                               ///< gives the $\theta$ volume = fabs(cos(th_m) - cos(th_p))

  std::vector<real> xbeg;             ///< Beginning of active domain in datablock
  std::vector<real> xend;             ///< End of active domain in datablock

  IdefixArray3D<real> dV;                ///< cell volume
  std::vector<IdefixArray3D<real>> A;    ///< cell right interface area

  std::vector<int> np_tot;        ///< total number of grid points in datablock
  std::vector<int> np_int;        ///< active number of grid points in datablock (excl. ghost cells)

  std::vector<int> nghost;               ///< number of ghost cells at each boundary
  std::vector<BoundaryType> lbound;      ///< Boundary condition to the left
  std::vector<BoundaryType> rbound;      ///< Boundary condition to the right

  bool haveAxis{false};       ///< DataBlock contains points on the axis and a special treatment
                              ///< has been required for these.

  std::vector<int> beg;       ///< First local index of the active domain
  std::vector<int> end;       ///< Last local index of the active domain+1

  std::vector<int> gbeg;      ///< First global index of the active domain of this datablock
  std::vector<int> gend;      ///< Last global index of the active domain of this datablock

  real dt;                     ///< Current timestep
  real t;                      ///< Current time

  Grid *mygrid;                ///< Parent grid object

  StateContainer states;       ///< conservative state of the datablock (contains references to dedicated objects)

  Hydro hydro;                  ///< The Hydro object attached to this datablock

  void InitFromGrid(Grid &, Input &); ///< init from a Grid object
  void MakeGeometry();                ///< Compute geometrical terms
  void DumpToFile(std::string);   ///< Dump current datablock to a file for inspection
  int CheckNan();                 ///< Return the number of cells which have Nans

  bool rklCycle{false};           ///<  // Set to true when we're inside a RKL call

  void EvolveStage();             ///< Evolve this DataBlock by dt
  void SetBoundaries();       ///< Enforce boundary conditions to this datablock
  void ShowConfig();              ///< Show the datablock's configuration

  void ResetStage();              ///< Reset the variables needed at each major integration Stage

  DataBlock();

  // Do we use fargo-like scheme ? (orbital advection)
  bool haveFargo{false};
  Fargo fargo;

  // Do we have Gravity ?
  bool haveGravity{false};
  Gravity gravity;

 private:
  void WriteVariable(FILE* , int , int *, char *, void*);

  template<int dir> void LoopDir();     ///< // recursive loop on dimensions
};

#endif // DATABLOCK_DATABLOCK_HPP_
