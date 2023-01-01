// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef FLUID_VISCOSITY_HPP_
#define FLUID_VISCOSITY_HPP_

#include "idefix.hpp"
#include "input.hpp"
#include "grid.hpp"
#include "fluid_defs.hpp"



// Forward class hydro declaration
template <typename Phys> class Fluid;
class DataBlock;

using ViscousDiffusivityFunc = void (*) (DataBlock &, const real t,
                                         IdefixArray3D<real> &, IdefixArray3D<real> &);

class Viscosity {
 public:
  template <typename Phys>
  Viscosity(Input &, Grid &, Fluid<Phys> *); ;  // constructor

  void ShowConfig();                    // print configuration

  void AddViscousFlux(int, const real, const IdefixArray4D<real> &);

  // Enroll user-defined viscous diffusivity
  void EnrollViscousDiffusivity(ViscousDiffusivityFunc);

  IdefixArray4D<real> viscSrc;  // Source terms of the viscous operator
  IdefixArray3D<real> eta1Arr;
  IdefixArray3D<real> eta2Arr;

  // pre-computed geometrical factors in non-cartesian geometry
  IdefixArray1D<real> one_dmu;

 private:
  DataBlock* data;

  ParabolicModuleStatus status;

  // type of viscosity function
  HydroModuleStatus haveViscosity{Disabled};
  ViscousDiffusivityFunc viscousDiffusivityFunc;

  IdefixArray4D<real> &Vc;
  IdefixArray3D<real> &dMax;

  // constant diffusion coefficient (when needed)
  real eta1, eta2;
};

#include "fluid.hpp"

template<typename Phys>
Viscosity::Viscosity(Input &input, Grid &grid, Fluid<Phys> *hydroin): 
                      Vc{hydroin->Vc}, 
                      dMax{hydroin->dMax},
                      status{hydroin->viscosityStatus} {

  idfx::pushRegion("Viscosity::Viscosity");
  // Save the parent hydro object
  this->data = hydroin->data;

  if(input.CheckEntry("Hydro","viscosity")>=0) {
    if(input.Get<std::string>("Hydro","viscosity",1).compare("constant") == 0) {
        this->eta1 = input.Get<real>("Hydro","viscosity",2);
        // second viscosity?
        this->eta2 = input.GetOrSet<real>("Hydro","viscosity",3, 0.0);
        this->haveViscosity = Constant;
      } else if(input.Get<std::string>("Hydro","viscosity",1).compare("userdef") == 0) {
        this->haveViscosity = UserDefFunction;
        this->eta1Arr = IdefixArray3D<real>("ViscosityEta1Array",data->np_tot[KDIR],
                                                                 data->np_tot[JDIR],
                                                                 data->np_tot[IDIR]);
        this->eta2Arr = IdefixArray3D<real>("ViscosityEta1Array",data->np_tot[KDIR],
                                                                 data->np_tot[JDIR],
                                                                 data->np_tot[IDIR]);

      } else {
        IDEFIX_ERROR("Unknown viscosity definition in idefix.ini. "
                     "Can only be constant or userdef.");
      }
  } else {
    IDEFIX_ERROR("I cannot create a Viscosity object without viscosity defined"
                   "in the .ini file");
  }

  // Allocate and fill arrays when needed
  #if GEOMETRY != CARTESIAN
    one_dmu = IdefixArray1D<real>("Viscosity_1dmu", data->np_tot[JDIR]);
    IdefixArray1D<real> dmu = one_dmu;
    IdefixArray1D<real> th = data->x[JDIR];
    idefix_for("ViscousInitGeometry",1,data->np_tot[JDIR],
      KOKKOS_LAMBDA(int j) {
        real scrch =  FABS((1.0-cos(th(j)))*(sin(th(j)) >= 0.0 ? 1.0:-1.0)
                     -(1.0-cos(th(j-1))) * (sin(th(j-1)) > 0.0 ? 1.0:-1.0));
        dmu(j) = 1.0/scrch;
      });
  #endif
  viscSrc = IdefixArray4D<real>("Viscosity_source", COMPONENTS, data->np_tot[KDIR],
                                                                data->np_tot[JDIR],
                                                                data->np_tot[IDIR]);
  idfx::popRegion();
}
#endif // FLUID_VISCOSITY_HPP_
