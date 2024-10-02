// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#ifndef FLUID_EVOLVESTAGE_HPP_
#define FLUID_EVOLVESTAGE_HPP_

#include "fluid.hpp"
#include "riemannSolver.hpp"
#include "forcing.hpp"
template<typename Phys>
template<int dir>
void Fluid<Phys>::LoopDir(const real t, const real dt) {
    // Step 2: compute the intercell flux with our Riemann solver, store the resulting InvDt
    this->rSolver->template CalcFlux<dir>(this->FluxRiemann);

    // Step 2.5: compute intercell parabolic flux when needed
    if(haveExplicitParabolicTerms) CalcParabolicFlux<dir>(t);

    // If we have tracers, compute the tracer intercell flux
    if(haveTracer) {
      this->tracer->template CalcFlux<dir, Phys>(this->FluxRiemann);
    }

    // Step 3: compute the resulting evolution of the conserved variables, stored in Uc
    CalcRightHandSide<dir>(t,dt);
    if(haveTracer) {
      this->tracer->template CalcRightHandSide<dir, Phys>(this->FluxRiemann,t ,dt);
    }

    // Recursive: do next dimension
    if constexpr (dir+1 < DIMENSIONS) LoopDir<dir+1>(t, dt);
}



// Evolve one step forward in time of hydro
template<typename Phys>
void Fluid<Phys>::EvolveStage(const real t, const real dt) {
  idfx::pushRegion("Fluid::EvolveStage");
  // Compute current when needed
  if(needExplicitCurrent) CalcCurrent();

  if(hallStatus.status == UserDefFunction) {
    if(hallDiffusivityFunc)
      hallDiffusivityFunc(*data, t, xHall);
    else
      IDEFIX_ERROR("No user-defined Hall diffusivity function has been enrolled");
  }

  if constexpr(Phys::eos) {
    eos->Refresh(*data, t);
  }

  // Loop on all of the directions
  LoopDir<IDIR>(t,dt);

  // Step 4: add source terms to the conserved variables (curvature, rotation, etc)
  if(haveSourceTerms) AddSourceTerms(t, dt);

  // Step 5: add drag when needed
  if(haveDrag) drag->AddDragForce(dt);

  if constexpr(Phys::mhd) {
    #if DIMENSIONS >= 2
      // Compute the field evolution according to CT
      emf->CalcCornerEMF(t);
      if(resistivityStatus.isExplicit || ambipolarStatus.isExplicit) {
        emf->CalcNonidealEMF(t);
      }
      emf->EnforceEMFBoundary();
      #ifdef EVOLVE_VECTOR_POTENTIAL
        emf->EvolveVectorPotential(dt, Ve);
        emf->ComputeMagFieldFromA(Ve, Vs);
      #else
        emf->EvolveMagField(t, dt, Vs);
      #endif

      boundary->ReconstructVcField(Uc);
    #endif
  }

  idfx::popRegion();
}

// Evolve one step forward in time of forcing
template<typename Phys>
void Fluid<Phys>::EvolveForcing(const real t, const real dt) {
  idfx::pushRegion("Fluid::EvolveForcing");

  data->forcing->ComputeForcing(dt);

  // Loop on all of the directions
  LoopForcingDir<IDIR>(t,dt);

  idfx::popRegion();
}

template<typename Phys>
template<int dir>
void Fluid<Phys>::LoopForcingDir(const real t, const real dt) {
//  IdefixArray3D<real> dV;
//  IdefixArray4D<real> Flux;
//  IdefixArray4D<real> forcingTerm;
//
//  IdefixArray1D<real> rt;
//  IdefixArray1D<real> dmu;
//  IdefixArray1D<real> dx;
//  IdefixArray1D<real> dx2;
//  real rhs[Phys::nvar]; //WARNING is it zero uniformly when called like that?
//  forcingTerm = data->forcing->forcingTerm;

  IdefixArray3D<real> dV = data->dV;
  IdefixArray4D<real> Flux = FluxRiemann;
  IdefixArray4D<real> forcingTerm = data->forcing->forcingTerm;

  IdefixArray1D<real> rt = data->rt;
  IdefixArray1D<real> dmu = data->dmu;
  IdefixArray1D<real> dx = data->dx[dir];
  IdefixArray1D<real> dx2 = data->dx[JDIR];
  IdefixArray1D<real> rhs = IdefixArray1D<real>("RhsForcing", data->np_tot[dir]);

//  KOKKOS_INLINE_FUNCTION void operator() (const int k, const int j,  const int i) const {
  idefix_for("ComputeForcing", 0, data->np_tot[KDIR], 0, data->np_tot[JDIR], 0, data->np_tot[IDIR],                KOKKOS_LAMBDA (int k, int j, int i) {

    const int ioffset = (dir==IDIR) ? 1 : 0;
    const int joffset = (dir==JDIR) ? 1 : 0;
    const int koffset = (dir==KDIR) ? 1 : 0;

    real dtdV=dt / dV(k,j,i);

    // elmentary length for gradient computations
    const int ig = ioffset*i + joffset*j + koffset*k;
    real dl = dx(ig);
    #if GEOMETRY == POLAR
      if constexpr (dir==JDIR)
        dl = dl*x1(i);

    #elif GEOMETRY == SPHERICAL
      if constexpr(dir==JDIR)
        dl = dl*rt(i);
      else if constexpr(dir==KDIR)
          dl = dl*rt(i)*dmu(j)/dx2(j);
    #endif


      rhs(MX1+dir) += dt * Vc(RHO,k,j,i) * forcingTerm(dir,k,j,i);
      if constexpr(Phys::pressure) {
        //  rho * v . f, where rhov is taken as a  volume average of Flux(RHO)
        rhs(ENG) += HALF_F * dtdV * dl *
                      (Flux(RHO,k,j,i) + Flux(RHO, k+koffset, j+joffset, i+ioffset)) *
                        forcingTerm(dir,k,j,i);
      } // Pressure

      // Particular cases if we do not sweep all of the components
      #if DIMENSIONS == 1 && COMPONENTS > 1
        EXPAND(                                                           ,
                  rhs(MX2) += dt * Vc(RHO,k,j,i) * forcingTerm(JDIR,k,j,i);   ,
                  rhs(MX3) += dt * Vc(RHO,k,j,i) * forcingTerm(KDIR,k,j,i);    )
        if constexpr(Phys::pressure) {
          rhs(ENG) += dt * (EXPAND( ZERO_F                                              ,
                                    + Vc(RHO,k,j,i) * Vc(VX2,k,j,i) * forcingTerm(JDIR,k,j,i)   ,
                                    + Vc(RHO,k,j,i) * Vc(VX3,k,j,i) * forcingTerm(KDIR,k,j,i) ));
        }
      #endif
      #if DIMENSIONS == 2 && COMPONENTS == 3
        // Only add this term once!
        if constexpr (dir==JDIR) {
          rhs(MX3) += dt * Vc(RHO,k,j,i) * forcingTerm(KDIR,k,j,i);
          if constexpr(Phys::pressure) {
            rhs(ENG) += dt * Vc(RHO,k,j,i) * Vc(VX3,k,j,i) * forcingTerm(KDIR,k,j,i);
          }
        }
      #endif

  #pragma unroll
  for(int nv = 0 ; nv < Phys::nvar ; nv++) {
    Uc(nv,k,j,i) = Uc(nv,k,j,i) + rhs(nv);
  }

  // WARNING DO SOMETHING FOR THE TIMESTEP DT?

  });

  // Recursive: do next dimension
  if constexpr (dir+1 < DIMENSIONS) LoopDir<dir+1>(t, dt);
}
#endif //FLUID_EVOLVESTAGE_HPP_
