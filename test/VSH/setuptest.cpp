#include "idefix.hpp"
#include "setup.hpp"
#include <vector>
#include <fstream>
#include <random>
#include <cmath>


// Initialisation routine. Can be used to allocate Arrays or variables which are used later on
Setup::Setup(Input &input, Grid &grid, DataBlock &data, Output &output) {
}

// This routine initializes the flow
// Note that data is on the device.
// One can therefore define locally
// a datahost and sync it, if needed
void Setup::InitFlow(DataBlock &data) {
  // Create a host copy
//  DataBlockHost d(data);
  DataBlockHost d(data, 1);
 
  // Shallow copy the global array
  for(int k = 0; k < d.np_tot[KDIR] ; k++) {
    for(int j = 0; j < d.np_tot[JDIR] ; j++) {
      for(int i = 0; i < d.np_tot[IDIR] ; i++) {

        real x1=d.x[IDIR](i);
        real x2=d.x[JDIR](j);
        real x3=d.x[KDIR](k);
        d.Vc(RHO,k,j,i) = 1.;
        d.Vc(PRS,k,j,i) = 1.;

//        d.Vc(VX1,k,j,i) = ZERO_F;
//          d.Vc(VX2,k,j,i) = 1.;
//        d.Vc(VX2,k,j,i) = ZERO_F;
//        d.Vc(VX3,k,j,i) = ZERO_F;
//        d.Vc(VX1,k,j,i) = vsh.Ylm_r(1,0,k,j);
//        d.Vc(VX1,k,j,i) = vsh.Ylm_r(1,0,k,j)/pow(x1,2.);
//        d.Vc(VX1,k,j,i) = amp*d.vsh->Ylm_r(1,0,k,j);
//        d.Vc(VX2,k,j,i) = amp*d.vsh->Slm_th(1,1,k,j);
//        d.Vc(VX3,k,j,i) = amp*d.vsh->Slm_phi(1,1,k,j);
//        d.Vc(VX1,k,j,i) = amp*;
//        d.Vc(VX2,k,j,i) = amp*d.vsh->Slm_th(1,1,k,j);
//        d.Vc(VX3,k,j,i) = amp*d.vsh->Slm_phi(1,1,k,j);
//        d.Vc(VX1,k,j,i) = amp*normal_distrib(generator);
//        d.Vc(VX2,k,j,i) = amp*normal_distrib(generator);
//        d.Vc(VX3,k,j,i) = amp*normal_distrib(generator);



















d.Vc(VX1,k,j,i) = d.Ylm_r(3,3,k,j).real();
d.Vc(VX2,k,j,i) = d.Tlm_th(3,3,k,j).real();
d.Vc(VX3,k,j,i) = d.Tlm_phi(3,3,k,j).real();

//          d.Vs(BX1s,k,j,i) = ZERO_F;
//          d.Vs(BX2s,k,j,i) = amp*d.vsh->Tlm_ths(2,2,k,j);
//          d.Vs(BX3s,k,j,i) = amp*d.vsh->Tlm_phis(2,2,k,j); 
      }
    }
  } 
  d.SyncToDevice();
}

