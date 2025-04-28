#include "idefix.hpp"
#include "setup.hpp"
#include "analysis.hpp"
#include <vector>
#include <fstream>
#include <random>
//#include <cmath>

static real rminGlob;
static real rmaxGlob;
static real thetaminGlob;
static real thetamaxGlob;
static real phiminGlob;
static real phimaxGlob;
static int nrGlob;
static int nthGlob;
static int nphiGlob;

static real ampGlob;

enum DiffusionType {conductivity, diffusivity, spitzer};
DiffusionType diffTypeGlob;

enum MagType {radial, azimuthal};
MagType magTypeGlob; 

Analysis *analysis;
void AnalysisFunction(DataBlock &data) {
  analysis->PerformAnalysis(data);
}

// Initialisation routine. Can be used to allocate Arrays or variables which are used later on
Setup::Setup(Input &input, Grid &grid, DataBlock &data, Output &output) {
  std::string outputFolder = input.GetOrSet<std::string>("Output","folder",0,"output");
  if(idfx::prank==0) {
    if(!fs::is_directory(outputFolder)) {
      try {
        if(!fs::create_directory(outputFolder)) {
          std::stringstream msg;
          msg << "Cannot create directory " << outputFolder << std::endl;
          IDEFIX_ERROR(msg);
        }
      } catch(std::exception &e) {
        std::stringstream msg;
        msg << "Cannot create directory " << outputFolder << std::endl;
        msg << e.what();
        IDEFIX_ERROR(msg);
      }
    }
  }
  analysis = new Analysis(input, grid, data, output,std::string(outputFolder+"/timevol.dat"));
  output.EnrollAnalysis(&AnalysisFunction);
  // Reset analysis if required
  if(!input.restartRequested) {
    analysis->ResetAnalysis();
  }

  ampGlob = input.Get<real>("Setup","amp",0);

  rminGlob = input.Get<real>("Grid","X1-grid",1);
  rmaxGlob = input.Get<real>("Grid","X1-grid",4);
  thetaminGlob = input.Get<real>("Grid","X2-grid",1);
  thetamaxGlob = input.Get<real>("Grid","X2-grid",4);
  phiminGlob = input.Get<real>("Grid","X3-grid",1);
  phimaxGlob = input.Get<real>("Grid","X3-grid",4);
  nrGlob = input.Get<real>("Grid","X1-grid",2);
  nthGlob = input.Get<real>("Grid","X2-grid",2);
  nphiGlob = input.Get<real>("Grid","X3-grid",2);
  real rmin = rminGlob;
  real rmax = rmaxGlob;
}

// This routine initializes the flow
// Note that data is on the device.
// One can therefore define locally
// a datahost and sync it, if needed
void Setup::InitFlow(DataBlock &data) {
  // Create a host copy
  DataBlockHost d(data);
 
  real amp = ampGlob;
  real kr = M_PI/rminGlob;
//  std::default_random_engine generator(idfx::prank);
//  std::default_random_engine generatorB(idfx::prank);
  std::default_random_engine generator(0);
  std::default_random_engine generatorB(0);
  std::normal_distribution<real> normal_distrib(0.,1.);

  int kghost = d.nghost[KDIR];
  int jghost = d.nghost[JDIR];
  int ighost = d.nghost[IDIR];
  int nphitot = d.np_tot[KDIR];
  int nthtot = d.np_tot[JDIR];
  int nrtot = d.np_tot[IDIR];
  IdefixHostArray3D<real> Vx1Array = IdefixHostArray3D<real>("Vx1Array", nphitot, nthtot, nrtot);
  IdefixHostArray3D<real> Vx2Array = IdefixHostArray3D<real>("Vx2Array", nphitot, nthtot, nrtot);
  IdefixHostArray3D<real> Vx3Array = IdefixHostArray3D<real>("Vx3Array", nphitot, nthtot, nrtot);
  IdefixHostArray3D<real> Bx1Array = IdefixHostArray3D<real>("Bx1Array", nphitot, nthtot, nrtot);
  IdefixHostArray3D<real> Bx2Array = IdefixHostArray3D<real>("Bx2Array", nphitot, nthtot, nrtot);
  IdefixHostArray3D<real> Bx3Array = IdefixHostArray3D<real>("Bx3Array", nphitot, nthtot, nrtot);
  for(int i = 0; i < nrtot ; i++) {
    for(int j = 0; j < nthtot ; j++) {
      for(int k = 0; k < nphitot ; k++) {
        Vx1Array(k,j,i) = ZERO_F; 
        Vx2Array(k,j,i) = ZERO_F; 
        Vx3Array(k,j,i) = ZERO_F;
        Bx1Array(k,j,i) = ZERO_F;
        Bx2Array(k,j,i) = ZERO_F;
        Bx3Array(k,j,i) = ZERO_F;
      }
    }
  }
  int kbeg = d.gbeg[KDIR] - kghost;
  int jbeg = d.gbeg[JDIR] - jghost;
  int ibeg = d.gbeg[IDIR] - ighost;
  int kend = d.gend[KDIR];
  int jend = d.gend[JDIR];
  int iend = d.gend[IDIR];
  for(int i = 0; i < nrGlob+2*ighost ; i++) {
    for(int j = 0; j < nthGlob+2*jghost ; j++) {
      for(int k = 0; k < nphiGlob+2*kghost ; k++) {
        real current_Vx1 = amp*normal_distrib(generator);
        real current_Vx2 = amp*normal_distrib(generator);
        real current_Vx3 = amp*normal_distrib(generator);
        real current_Bx1 = amp*normal_distrib(generatorB);
        real current_Bx2 = amp*normal_distrib(generatorB);
        real current_Bx3 = amp*normal_distrib(generatorB);

        if (i >= ibeg and i < iend and j >= jbeg and j < jend and k >= kbeg and k < kend) {
          Vx1Array(k-kbeg,j-jbeg,i-ibeg) = current_Vx1;
          Vx2Array(k-kbeg,j-jbeg,i-ibeg) = current_Vx2;
          Vx3Array(k-kbeg,j-jbeg,i-ibeg) = current_Vx3;
          Bx1Array(k-kbeg,j-jbeg,i-ibeg) = current_Bx1;
          Bx2Array(k-kbeg,j-jbeg,i-ibeg) = current_Bx2;
          Bx3Array(k-kbeg,j-jbeg,i-ibeg) = current_Bx3;
        }
      }
    }
  }

  for(int i = 0; i < d.np_tot[IDIR] ; i++) {
    real x1=d.x[IDIR](i);
    for(int j = 0; j < d.np_tot[JDIR] ; j++) {
      for(int k = 0; k < d.np_tot[KDIR] ; k++) {
        d.Vc(RHO,k,j,i) = 1.;
//        d.Vc(PRS,k,j,i) = 1.;

        real x1=d.x[IDIR](i);
        real x1l=d.xl[IDIR](i);
        real x2=d.x[JDIR](j);
        real x2l=d.xl[JDIR](j);
        real x3=d.x[KDIR](k);
        real x3l=d.xl[KDIR](k);

        d.Vc(VX1,k,j,i) = Vx1Array(k,j,i);
        d.Vc(VX2,k,j,i) = Vx2Array(k,j,i);
        #if COMPONENTS == 3
          d.Vc(VX3,k,j,i) = Vx3Array(k,j,i);
        #endif //COMPONENTS == 3
      }
    }
  }
  d.SyncToDevice();
}

