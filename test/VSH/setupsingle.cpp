#include "idefix.hpp"
#include "setup.hpp"

// Initialisation routine. Can be used to allocate Arrays or variables which are used later on
Setup::Setup(Input &input, Grid &grid, DataBlock &data, Output &output) {
}

// This routine initializes the flow
// Note that data is on the device.
// One can therefore define locally
// a datahost and sync it, if needed
void Setup::InitFlow(DataBlock &data) {
    // Create a host copy
  DataBlockHost d(data, 1);
//  d.SyncToDevice();
//  std::default_random_engine generator;

  int countr = 0;
  int countth = 0;
  int countphi = 0;
  for(int k = 2; k < d.np_tot[KDIR]-2 ; k++) {
    for(int j = 2; j < d.np_tot[JDIR]-2 ; j++) {
      for(int i = 2; i < d.np_tot[IDIR]-2 ; i++) {
        real x1=d.x[IDIR](i);
        real x2=d.x[JDIR](j);
        real x3=d.x[KDIR](k);
        d.Vc(RHO,k,j,i) = 1.;
        d.Vc(PRS,k,j,i) = 1.;
        
//        d.Vc(VX1,k,j,i) = data.Ylm_r(2,2,k,j);
//        d.Vc(VX2,k,j,i) = data.Tlm_th(1,1,k,j);
//        d.Vc(VX3,k,j,i) = data.Tlm_phi(1,1,k,j);
        d.Vc(VX1,k,j,i) = d.Ylm_r(2,2,k,j).real();
        d.Vc(VX2,k,j,i) = d.Tlm_th(1,1,k,j).real();
        d.Vc(VX3,k,j,i) = d.Tlm_phi(1,1,k,j).real();
        if (FABS(d.Vc(VX3,k,j,i) - sqrt(2.*3./(8.*M_PI))*COS(x2)*COS(x3)) > 1e-6) {
          countphi ++;
        }
        if (FABS(d.Vc(VX2,k,j,i) - sqrt(2.*3./(8.*M_PI))*SIN(x3)) > 1e-6) {
          countth ++;
        }
        if (FABS(d.Vc(VX1,k,j,i) - 0.5*sqrt(15./(2.*M_PI))*pow(SIN(x2),2)*COS(2.*x3)) > 1e-6) {
          countr ++;
        }
      }
    }
  } 
  if (countphi==0 & countth==0 & countr==0) {
    std::cout << "SUCCESS: no discrepancy has been found." << std::endl;
  } else {
    std::cout << "FAILURE: countphi=" << countphi << ", countth=" << countth << " and countr=" << countr << std::endl;
    IDEFIX_ERROR("Error");
  }
  d.SyncToDevice();
}

