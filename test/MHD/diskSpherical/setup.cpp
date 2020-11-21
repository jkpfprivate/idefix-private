#include "idefix.hpp"
#include "setup.hpp"


real epsilonGlob;
real betaGlob;
real HidealGlob;
real AmMidGlob;
real gammaGlob;
real densityFloorGlob;

/*********************************************/
/**
Customized random number generator
Allow one to have consistant random numbers
generators on different architectures.
**/
/*********************************************/
real randm(void) {
    const int a    =    16807;
    const int m =    2147483647;
    static int in0 = 13763 + 2417*idfx::prank;
    int q;

    /* find random number  */
    q= (int) fmod((double) a * in0, m);
    in0=q;

    return((real) ((double) q/(double)m));
}




void MySourceTerm(DataBlock &data, const real t, const real dtin) {
  IdefixArray4D<real> Vc = data.Vc;
  IdefixArray4D<real> Uc = data.Uc;
  IdefixArray1D<real> x1=data.x[IDIR];
  IdefixArray1D<real> x2=data.x[JDIR];
  real epsilon = epsilonGlob;
  real tauGlob=0.1;
  real gamma_m1=gammaGlob-1.0;
  real dt=dtin;
  idefix_for("MySourceTerm",0,data.np_tot[KDIR],0,data.np_tot[JDIR],0,data.np_tot[IDIR],
              KOKKOS_LAMBDA (int k, int j, int i) {
                real r=x1(i);
                real th=x2(j);
                real z=r*cos(th);
                real R=r*sin(th);
                real cs2, tau;
                if(R>1.0) {
                    real csdisk = epsilon/sqrt(R);
                    cs2=(csdisk*csdisk);
                    tau= tauGlob*sqrt(R);
                }
                else {
                    real csdisk = epsilon;
                    cs2=(csdisk*csdisk);
                    tau=tauGlob;
                }
                // Cooling /heatig function
                real Ptarget = cs2*Vc(RHO,k,j,i);

                Uc(ENG,k,j,i) += -dt*(Vc(PRS,k,j,i)-Ptarget)/(tau*gamma_m1);

		// Velocity relaxation
		if(R<1.2) {
			Uc(MX1,k,j,i) += -dt*(Vc(VX1,k,j,i)*Vc(RHO,k,j,i));
			Uc(MX2,k,j,i) += -dt*(Vc(VX2,k,j,i)*Vc(RHO,k,j,i));
		}

});


}

void EmfBoundary(DataBlock& data, const real t) {
    IdefixArray3D<real> Ex1 = data.emf.ex;
    IdefixArray3D<real> Ex2 = data.emf.ey;
    IdefixArray3D<real> Ex3 = data.emf.ez;
    if(data.lbound[IDIR] == userdef) {

        int ighost = data.nghost[IDIR];

        idefix_for("EMFBoundary",0,data.np_tot[KDIR],0,data.np_tot[JDIR],0,ighost+1,
                    KOKKOS_LAMBDA (int k, int j, int i) {
            Ex3(k,j,i) = ZERO_F;
        });
    }
    if(data.lbound[JDIR] == userdef) {
        int jghost = data.nghost[JDIR];
        //printf("I'mbeing called\n");
        idefix_for("EMFBoundary",0,data.np_tot[KDIR],0,data.np_tot[IDIR],
                    KOKKOS_LAMBDA (int k, int i) {
            Ex3(k,jghost,i) = ZERO_F;
        });
    }
    if(data.rbound[JDIR] == userdef) {
        int jghost = data.end[JDIR];
        //printf("I'mbeing called\n");
        idefix_for("EMFBoundary",0,data.np_tot[KDIR],0,data.np_tot[IDIR],
                    KOKKOS_LAMBDA (int k, int i) {
            Ex3(k,jghost,i) = ZERO_F;
        });
    }
}

void InternalBoundary(DataBlock& data, const real t) {
  IdefixArray4D<real> Vc = data.Vc;
  IdefixArray4D<real> Vs = data.Vs;
  IdefixArray1D<real> x1=data.x[IDIR];
  IdefixArray1D<real> x2=data.x[JDIR];

  real vAmax=10.0;
  real densityFloor = densityFloorGlob;
  idefix_for("InternalBoundary",0,data.np_tot[KDIR],0,data.np_tot[JDIR],0,data.np_tot[IDIR],
              KOKKOS_LAMBDA (int k, int j, int i) {
                real b2=EXPAND(Vc(BX1,k,j,i)*Vc(BX1,k,j,i) , +Vc(BX2,k,j,i)*Vc(BX2,k,j,i), +Vc(BX3,k,j,i)*Vc(BX3,k,j,i) ) ;
                real va2=b2/Vc(RHO,k,j,i);
                real myMax=vAmax;
                //if(x1(i)<1.1) myMax=myMax/50.0;
                if(va2>myMax*myMax) {
                  real T = Vc(PRS,k,j,i)/Vc(RHO,k,j,i);
                  Vc(RHO,k,j,i) = b2/(myMax*myMax);
                  Vc(PRS,k,j,i) = T*Vc(RHO,k,j,i);
                }
                if(Vc(RHO,k,j,i) < densityFloor) {
                  real T= Vc(PRS,k,j,i)/Vc(RHO,k,j,i);
                  Vc(RHO,k,j,i)=densityFloor;
                  Vc(PRS,k,j,i)=T*Vc(RHO,k,j,i);
                }
                /*
                  real R = x1(i)*sin(x2(j));
                  if(R<1.0) {
                      Vc(VX1,k,j,i) = ZERO_F;
                      Vc(VX2,k,j,i) = ZERO_F;
                      Vc(VX3,k,j,i) = R;
                  }*/
              });

}
// User-defined boundaries
void UserdefBoundary(DataBlock& data, int dir, BoundarySide side, real t) {

    if( (dir==IDIR) && (side == left)) {
        IdefixArray4D<real> Vc = data.Vc;
        IdefixArray4D<real> Vs = data.Vs;
        IdefixArray1D<real> x1 = data.x[IDIR];
        IdefixArray1D<real> x2 = data.x[JDIR];

        int ighost = data.nghost[IDIR];
        real Omega=1.0;
        idefix_for("UserDefBoundary",0,data.np_tot[KDIR],0,data.np_tot[JDIR],0,ighost,
                    KOKKOS_LAMBDA (int k, int j, int i) {
                        real R=x1(i)*sin(x2(j));
                        real z=x1(i)*cos(x2(j));

                        Vc(RHO,k,j,i) = Vc(RHO,k,j,ighost);
                        Vc(PRS,k,j,i) = Vc(PRS,k,j,ighost);
                        if(Vc(VX1,k,j,ighost)>=ZERO_F) Vc(VX1,k,j,i) = - Vc(VX1,k,j,2*ighost-i);
			                   else Vc(VX1,k,j,i) = Vc(VX1,k,j,ighost);
                        Vc(VX2,k,j,i) = Vc(VX2,k,j,ighost);
                        Vc(VX3,k,j,i) = R*Omega;
                        Vs(BX2s,k,j,i) = Vs(BX2s,k,j,ighost);
                        Vc(BX3,k,j,i) = ZERO_F;

                    });

    }

    if( dir==JDIR) {
        IdefixArray4D<real> Vc = data.Vc;
        IdefixArray4D<real> Vs = data.Vs;
        int jghost;
        int jbeg,jend;
        if(side == left) {
            jghost = data.beg[JDIR];
            jbeg = 0;
            jend = data.beg[JDIR];
            //return;
        }
        else if(side==right) {
            jghost = data.end[JDIR]-1;
            jbeg=data.end[JDIR];
            jend=data.np_tot[JDIR];
        }


        idefix_for("UserDefBoundary",0,data.np_tot[KDIR],jbeg,jend,0,data.np_tot[IDIR],
                    KOKKOS_LAMBDA (int k, int j, int i) {
                        Vc(RHO,k,j,i) = Vc(RHO,k,jghost,i);
                        Vc(PRS,k,j,i) = Vc(PRS,k,jghost,i);
                        Vc(VX1,k,j,i) = Vc(VX1,k,jghost,i);
                        Vc(VX2,k,j,i) = - Vc(VX2,k,2*jghost-j,i);
                        Vc(VX3,k,j,i) = ZERO_F;
                        Vs(BX1s,k,j,i) = Vs(BX1s,k,jghost,i);
                        Vc(BX3,k,j,i) = ZERO_F;

                    });


    }

}


void Potential(DataBlock& data, const real t, IdefixArray1D<real>& x1, IdefixArray1D<real>& x2, IdefixArray1D<real>& x3, IdefixArray3D<real>& phi) {

    idefix_for("Potential",0,data.np_tot[KDIR], 0, data.np_tot[JDIR], 0, data.np_tot[IDIR],
        KOKKOS_LAMBDA (int k, int j, int i) {
        phi(k,j,i) = -1.0/x1(i);
    });

}

// Default constructor
Setup::Setup() {}

// Initialisation routine. Can be used to allocate
// Arrays or variables which are used later on
Setup::Setup(Input &input, Grid &grid, DataBlock &data, Hydro &hydro) {
  // Set the function for userdefboundary
  hydro.EnrollUserDefBoundary(&UserdefBoundary);
  hydro.EnrollGravPotential(&Potential);
  hydro.EnrollUserSourceTerm(&MySourceTerm);
  hydro.EnrollInternalBoundary(&InternalBoundary);
  hydro.EnrollEmfBoundary(&EmfBoundary);
  gammaGlob=hydro.GetGamma();
  epsilonGlob = input.GetReal("Setup","epsilon",0);
  densityFloorGlob = input.GetReal("Setup","densityFloor",0);
}

// This routine initialize the flow
// Note that data is on the device.
// One can therefore define locally
// a datahost and sync it, if needed
void Setup::InitFlow(DataBlock &data) {
    // Create a host copy
    DataBlockHost d(data);

    // Make vector potential
    IdefixHostArray4D<real> A = IdefixHostArray4D<real>("Setup_VectorPotential", 3, data.np_tot[KDIR], data.np_tot[JDIR], data.np_tot[IDIR]);


    real x,y,z;

    real vphi,f,r,th;

    real epsilon=0.1;
    real beta=1000;

    for(int k = 0; k < d.np_tot[KDIR] ; k++) {
        for(int j = 0; j < d.np_tot[JDIR] ; j++) {
            for(int i = 0; i < d.np_tot[IDIR] ; i++) {
                r=d.x[IDIR](i);
                th=d.x[JDIR](j);
                real R=r*sin(th);
                real z=r*cos(th);
                real Vk=1.0/pow(R,0.5);

                real cs2=(epsilon*Vk)*(epsilon*Vk);

                d.Vc(RHO,k,j,i) = 1.0/(R*sqrt(R))*exp(1.0/(cs2)*(1/r-1/R));
                d.Vc(PRS,k,j,i) = cs2*d.Vc(RHO,k,j,i);
                d.Vc(VX1,k,j,i) = 0.0;
                d.Vc(VX2,k,j,i) = 1e-1*(0.5-randm());
                d.Vc(VX3,k,j,i) = Vk*sqrt(R/r-2.5*cs2);

                d.Vs(BX1s,k,j,i) = 0.0;
                d.Vs(BX2s,k,j,i) = 0.0;
                d.Vs(BX3s,k,j,i) = 0.0;

                real B0 = sqrt(2*cs2/(R*sqrt(R))/sqrt(beta));

                d.Vs(BX3s,k,j,i) = B0*cos(R/epsilon)*fmax(1-(z*z)/(4*R*R*epsilon*epsilon),ZERO_F);
                d.Vs(BX3s,k,j,i) *= fmax(tanh(10*(R-1.5)),ZERO_F);
                A(IDIR,k,j,i) = 0.0;
                A(JDIR,k,j,i) = 0.0;
                A(KDIR,k,j,i) = B0*epsilon*cos(R/epsilon)*fmax(1-(z*z)/(4*R*R*epsilon*epsilon),ZERO_F);
            }
        }
    }

    // Make the field from the vector potential
    //d.MakeVsFromAmag(A);

    // Send it all, if needed
    d.SyncToDevice();
}

// Analyse data to produce an output
void Setup::MakeAnalysis(DataBlock & data, real t) {

}




// Do a specifically designed user step in the middle of the integration
void ComputeUserStep(DataBlock &data, real t, real dt) {

}
