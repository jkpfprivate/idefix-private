// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2022 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2021 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

#include "vsh.hpp"
#include "dataBlockHost.hpp"

Vsh::Vsh(DataBlockHost *datain, int write) {
  this->nphi = datain->nphi_tot;
  this->nphi_proc = datain->np_int[KDIR];
  this->nphi_shtns = 2*nphi; //so to have the points at the right coordinates, because the phi-grid is equally spaced between 0 (included) and 2*M_PI (excluded)
  this->ntheta = datain->nth_tot;
  this->ntheta_proc = datain->np_int[JDIR];
  this->ntheta_shtns = 2*ntheta + 1;
  this->nr_proc = datain->np_int[IDIR];
  this->koffset =  datain->gbeg[KDIR] - datain->nghost[KDIR];
  this->joffset =  datain->gbeg[JDIR] - datain->nghost[JDIR];
  this->kghost = datain->nghost[KDIR];
  this->jghost = datain->nghost[JDIR];
  this->ighost = datain->nghost[IDIR];
  this->x1 = datain->x[IDIR];
  this->x1l = datain->xl[IDIR];

  this->lmax = datain->lmax;
  this->mmax = datain->mmax;

  this->write = write;

  this->jl = IdefixArray2D<real>::HostMirror ("jl", lmax, nr_proc+2*ighost);
  this->jls = IdefixArray2D<real>::HostMirror ("jl", lmax, nr_proc+2*ighost);
  this->Ylm_r = IdefixArray4D<real>::HostMirror ("Ylm_r", lmax, mmax, nphi_proc+2*kghost, ntheta_proc+2*jghost);
  this->Slm_th = IdefixArray4D<real>::HostMirror ("Slm_th", lmax, mmax, nphi_proc+2*kghost, ntheta_proc+2*jghost);
  this->Slm_phi = IdefixArray4D<real>::HostMirror ("Slm_phi", lmax, mmax, nphi_proc+2*kghost, ntheta_proc+2*jghost);
  this->Tlm_th = IdefixArray4D<real>::HostMirror ("Tlm_th", lmax, mmax, nphi_proc+2*kghost, ntheta_proc+2*jghost);
  this->Tlm_phi = IdefixArray4D<real>::HostMirror ("Tlm_phi", lmax, mmax, nphi_proc+2*kghost, ntheta_proc+2*jghost);

  this->Slm_ths = IdefixArray4D<real>::HostMirror ("Slm_ths", lmax, mmax, nphi_proc+2*kghost, ntheta_proc+2*jghost);
  this->Slm_phis = IdefixArray4D<real>::HostMirror ("Slm_phis", lmax, mmax, nphi_proc+2*kghost, ntheta_proc+2*jghost);
  this->Tlm_ths = IdefixArray4D<real>::HostMirror ("Tlm_ths", lmax, mmax, nphi_proc+2*kghost, ntheta_proc+2*jghost);
  this->Tlm_phis = IdefixArray4D<real>::HostMirror ("Tlm_phis", lmax, mmax, nphi_proc+2*kghost, ntheta_proc+2*jghost);

  for(int l = 0; l < this->lmax; l++) {
    for(int m = 0; m < this->mmax; m++) {
      for(int k = 0; k < nphi_proc+2*kghost; k++) {
        for(int j = 0; j < ntheta_proc+2*jghost; j++) {
          this->Ylm_r(l,m,k,j) = 0.;
          this->Slm_th(l,m,k,j) = 0.;
          this->Slm_phi(l,m,k,j) = 0.;
          this->Tlm_th(l,m,k,j) = 0.;
          this->Tlm_phi(l,m,k,j) = 0.;
        }
      }
    }
  }
} 

void Vsh::ShowConfig() {
  idfx::cout << "Vsh: ENABLED." << std::endl;
  idfx::cout << "Vsh: lmax=" << this->lmax << " and mmax=" << this->mmax << "." << std::endl;
}

void Vsh::Generatejl() {
  for(int l = 0; l < this->lmax; l++) {
    for(int i = 0; i < nr_proc+2*ighost; i++) {
      this->jl(l,i) = std::sph_bessel(l, x1(i));
      this->jls(l,i) = std::sph_bessel(l, x1l(i));
    }
  }
}

void Vsh::GenerateCellVsh() {
  shtns_cfg shtns;                // handle to a sht transform configuration
  int NLM;
  std::complex<real> *Ylm, *Slm, *Tlm;      // spherical harmonics coefficients (l,m space): complex numbers.
  real *Yr;                // real space : r
  real *Stheta, *Sphi;                // real space : theta,phi
  real *Ttheta, *Tphi;                // real space : theta,phi

  const int mres = 1;             // periodicity in phi (1 for full-sphere, 2 for half the sphere, 3 for 1/3, etc...)

  int verbose = (idfx::prank == 0) ? 1 : 0;
  shtns_verbose(verbose);                       // displays informations during initialization.
  shtns_use_threads(0);           // enable multi-threaded transforms (if supported).
//  shtns = shtns_init( sht_gauss, lmax, mmax, mres, ntheta, nphi );
  shtns = shtns_init( sht_reg_fast, lmax, mmax, mres, ntheta, nphi_shtns );
  //**** Default normalisation is SHTNS orthonormalized i.e. C_m Y_lm_ortho with C_m = 2 if m!=0 and 1 otherwise.
  //**** Default norm of SHTNS S/T_lm is sqrt(l(l_1)), so we divide by this quantity.

  NLM = shtns->nlm;
  Yr = (real *) shtns_malloc( NSPAT_ALLOC(shtns) * sizeof(real));
  Tphi = (real *) shtns_malloc( NSPAT_ALLOC(shtns) * sizeof(real));
  Ttheta = (real *) shtns_malloc( NSPAT_ALLOC(shtns) * sizeof(real));
  Sphi = (real *) shtns_malloc( NSPAT_ALLOC(shtns) * sizeof(real));
  Stheta = (real *) shtns_malloc( NSPAT_ALLOC(shtns) * sizeof(real));
  Ylm = (std::complex<real> *) shtns_malloc( NLM * sizeof(std::complex<real>));
  Slm = (std::complex<real> *) shtns_malloc( NLM * sizeof(std::complex<real>));
  Tlm = (std::complex<real> *) shtns_malloc( NLM * sizeof(std::complex<real>));

  for(int l = 0; l < this->lmax; l++) {
    int sqrt_ll = (l != 0) ? sqrt(l*(l+1)) : 1;
    for(int m = 0; m < std::min(this->mmax, l+1); m++) {
      LM_LOOP(shtns,  Ylm[lm]=0.0; Slm[lm]=0.0;  Tlm[lm] = 0.0; )
      Ylm[LM(shtns,l,m)] = 1.0;
      Slm[LM(shtns,l,m)] = 1.0/sqrt_ll;
      Tlm[LM(shtns,l,m)] = 1.0/sqrt_ll;
      SH_to_spat(shtns,Ylm,Yr);
      SHsph_to_spat(shtns,Slm,Stheta,Sphi);
      SHtor_to_spat(shtns,Tlm,Ttheta,Tphi);
      if (write==1 and idfx::prank==0) {
        std::string path;
        path = "datavsh/Yr"+std::to_string(l)+std::to_string(m);
        std::cout << path << std::endl;
        this->write_mx(path,Yr,nphi_shtns,ntheta);
        path = "datavsh/Stheta"+std::to_string(l)+std::to_string(m);
        std::cout << path << std::endl;
        this->write_mx(path,Stheta,nphi_shtns,ntheta);
        path = "datavsh/Sphi"+std::to_string(l)+std::to_string(m);
        std::cout << path << std::endl;
        this->write_mx(path,Sphi,nphi_shtns,ntheta);
        path = "datavsh/Ttheta"+std::to_string(l)+std::to_string(m);
        std::cout << path << std::endl;
        this->write_mx(path,Ttheta,nphi_shtns,ntheta);
        path = "datavsh/Tphi"+std::to_string(l)+std::to_string(m);
        std::cout << path << std::endl;
        this->write_mx(path,Tphi,nphi_shtns,ntheta);
      }

      int kbeg = (koffset==0) ? kghost : 0;
      int jbeg = (joffset==0) ? jghost : 0;
      int kend = (koffset+nphi_proc-kghost>=nphi) ? nphi_proc+kghost : nphi_proc+2*kghost;
      int jend = (joffset+ntheta_proc-jghost>=ntheta) ? ntheta_proc+jghost : ntheta_proc+2*jghost;
      for(int k = kbeg; k < kend; k++) {
        for(int j = jbeg; j < jend; j++) {
          this->Ylm_r(l,m,k,j) = Yr[(2*(k+koffset-kghost)+1)*ntheta+(j+joffset-jghost)];
          this->Slm_th(l,m,k,j) = Stheta[(2*(k+koffset-kghost)+1)*ntheta+(j+joffset-jghost)];
          this->Slm_phi(l,m,k,j) = Sphi[(2*(k+koffset-kghost)+1)*ntheta+(j+joffset-jghost)];
          this->Tlm_th(l,m,k,j) = Ttheta[(2*(k+koffset-kghost)+1)*ntheta+(j+joffset-jghost)];
          this->Tlm_phi(l,m,k,j) = Tphi[(2*(k+koffset-kghost)+1)*ntheta+(j+joffset-jghost)];
        }
      }
    }
  }
  shtns_free(Yr);
  shtns_free(Tphi);
  shtns_free(Ttheta);
  shtns_free(Sphi);
  shtns_free(Stheta);
  shtns_free(Ylm);
  shtns_free(Tlm);
  shtns_free(Slm);
  shtns_destroy(shtns);
}

//defined at the right interface between cells in the direction of the corresponding component, like the magnetic field BX1s,2s,3s
void Vsh::GenerateInterfaceVsh() { 
  shtns_cfg shtns;                // handle to a sht transform configuration
  int NLM;
  std::complex<real> *Slm, *Tlm;      // spherical harmonics coefficients (l,m space): complex numbers.
  real *Stheta, *Sphi;                // real space : theta,phi
  real *Ttheta, *Tphi;                // real space : theta,phi

  const int mres = 1;             // periodicity in phi (1 for full-sphere, 2 for half the sphere, 3 for 1/3, etc...)

  int verbose = (idfx::prank == 0) ? 1 : 0;
  shtns_verbose(verbose);                       // displays informations during initialization.
  shtns_use_threads(0);           // enable multi-threaded transforms (if supported).
//  shtns = shtns_init( sht_reg_fast, lmax, mmax, mres, ntheta_shtns, nphi_shtns );
  shtns = shtns_init( sht_reg_poles, lmax, mmax, mres, ntheta_shtns, nphi_shtns );

  NLM = shtns->nlm;
  Tphi = (real *) shtns_malloc( NSPAT_ALLOC(shtns) * sizeof(real));
  Ttheta = (real *) shtns_malloc( NSPAT_ALLOC(shtns) * sizeof(real));
  Sphi = (real *) shtns_malloc( NSPAT_ALLOC(shtns) * sizeof(real));
  Stheta = (real *) shtns_malloc( NSPAT_ALLOC(shtns) * sizeof(real));
  Slm = (std::complex<real> *) shtns_malloc( NLM * sizeof(std::complex<real>));
  Tlm = (std::complex<real> *) shtns_malloc( NLM * sizeof(std::complex<real>));

  for(int l = 0; l < this->lmax; l++) {
    int sqrt_ll = (l != 0) ? sqrt(l*(l+1)) : 1;
    for(int m = 0; m < this->mmax; m++) {
      LM_LOOP(shtns, Slm[lm]=0.0;  Tlm[lm] = 0.0; )
      Slm[LM(shtns,l,m)] = 1.0/sqrt_ll;
      Tlm[LM(shtns,l,m)] = 1.0/sqrt_ll;
      SHsph_to_spat(shtns,Slm,Stheta,Sphi);
      SHtor_to_spat(shtns,Tlm,Ttheta,Tphi);
      if (write==1 & idfx::prank==0) {
        std::string path;
        path = "datavsh/Sthetapoles"+std::to_string(l)+std::to_string(m);
        std::cout << path << std::endl;
        this->write_mx(path,Stheta,nphi_shtns,ntheta_shtns);
        path = "datavsh/Sphipoles"+std::to_string(l)+std::to_string(m);
        std::cout << path << std::endl;
        this->write_mx(path,Sphi,nphi_shtns,ntheta_shtns);
        path = "datavsh/Tthetapoles"+std::to_string(l)+std::to_string(m);
        std::cout << path << std::endl;
        this->write_mx(path,Ttheta,nphi_shtns,ntheta_shtns);
        path = "datavsh/Tphipoles"+std::to_string(l)+std::to_string(m);
        std::cout << path << std::endl;
        this->write_mx(path,Tphi,nphi_shtns,ntheta_shtns);
      }
      int kbeg = (koffset==0) ? kghost : 0;
      int kend = (koffset+nphi_proc>=nphi) ? nphi_proc+kghost : nphi_proc+2*kghost;
      int jbeg = (joffset==0) ? jghost : 0;
      int jend = (joffset+ntheta_proc>=ntheta) ? ntheta_proc+jghost : ntheta_proc+2*jghost;
      for(int k = kbeg; k < kend; k++) {
        for(int j = jbeg; j < jend; j++) {
          this->Slm_ths(l,m,k,j) = Stheta[(2*(k+koffset-kghost)+1)*ntheta_shtns+2*(j+joffset-jghost)];
          this->Slm_phis(l,m,k,j) = Sphi[(2*(k+koffset-kghost))*ntheta_shtns+2*(j+joffset-jghost)+1];
          this->Tlm_ths(l,m,k,j) = Ttheta[(2*(k+koffset-kghost)+1)*ntheta_shtns+2*(j+joffset-jghost)];
          this->Tlm_phis(l,m,k,j) = Tphi[(2*(k+koffset-kghost))*ntheta_shtns+2*(j+joffset-jghost)+1];
        }
      }
    }
  }
  shtns_free(Tphi);
  shtns_free(Ttheta);
  shtns_free(Sphi);
  shtns_free(Stheta);
  shtns_free(Tlm);
  shtns_free(Slm);
  shtns_destroy(shtns);
}

void Vsh::GenerateCellGhostVsh() {
  shtns_cfg shtns;                // handle to a sht transform configuration
  int NLM;
  std::complex<real> *Ylm, *Slm, *Tlm;      // spherical harmonics coefficients (l,m space): complex numbers.
  real *Yr;                // real space : r
  real *Stheta, *Sphi;                // real space : theta,phi
  real *Ttheta, *Tphi;                // real space : theta,phi

  const int mres = 1;             // periodicity in phi (1 for full-sphere, 2 for half the sphere, 3 for 1/3, etc...)

  int verbose = (idfx::prank == 0) ? 1 : 0;
  shtns_verbose(verbose);                       // displays informations during initialization.
  shtns_use_threads(0);           // enable multi-threaded transforms (if supported).
//  shtns = shtns_init( sht_gauss, lmax, mmax, mres, ntheta, nphi );
  shtns = shtns_init( sht_reg_fast, lmax, mmax, mres, ntheta, nphi_shtns );

  NLM = shtns->nlm;
  Yr = (real *) shtns_malloc( NSPAT_ALLOC(shtns) * sizeof(real));
  Tphi = (real *) shtns_malloc( NSPAT_ALLOC(shtns) * sizeof(real));
  Ttheta = (real *) shtns_malloc( NSPAT_ALLOC(shtns) * sizeof(real));
  Sphi = (real *) shtns_malloc( NSPAT_ALLOC(shtns) * sizeof(real));
  Stheta = (real *) shtns_malloc( NSPAT_ALLOC(shtns) * sizeof(real));
  Ylm = (std::complex<real> *) shtns_malloc( NLM * sizeof(std::complex<real>));
  Slm = (std::complex<real> *) shtns_malloc( NLM * sizeof(std::complex<real>));
  Tlm = (std::complex<real> *) shtns_malloc( NLM * sizeof(std::complex<real>));

  for(int l = 0; l < this->lmax; l++) {
    int sqrt_ll = (l != 0) ? sqrt(l*(l+1)) : 1;
    for(int m = 0; m < this->mmax; m++) {
      LM_LOOP(shtns,  Ylm[lm]=0.0; Slm[lm]=0.0;  Tlm[lm] = 0.0; )
      Ylm[LM(shtns,l,m)] = 1.0;
      Slm[LM(shtns,l,m)] = 1.0/sqrt_ll;
      Tlm[LM(shtns,l,m)] = 1.0/sqrt_ll;
      SH_to_spat(shtns,Ylm,Yr);
      SHsph_to_spat(shtns,Slm,Stheta,Sphi);
      SHtor_to_spat(shtns,Tlm,Ttheta,Tphi);

      if (koffset==0) {
        for(int k = 0; k < kghost; k++) {
          int jbeg = (joffset==0) ? jghost : 0;
          int jend = (joffset+ntheta_proc>=ntheta) ? ntheta_proc+jghost : ntheta_proc+2*jghost;
          for(int j = jbeg; j < jend; j++) {
            this->Ylm_r(l,m,k,j) = Yr[(2*(nphi-k+koffset-kghost)+1)*ntheta+(j+joffset-jghost)];
            this->Slm_th(l,m,k,j) = Stheta[(2*(nphi-k+koffset-kghost)+1)*ntheta+(j+joffset-jghost)];
            this->Slm_phi(l,m,k,j) = Sphi[(2*(nphi-k+koffset-kghost)+1)*ntheta+(j+joffset-jghost)];
            this->Tlm_th(l,m,k,j) = Ttheta[(2*(nphi-k+koffset-kghost)+1)*ntheta+(j+joffset-jghost)];
            this->Tlm_phi(l,m,k,j) = Tphi[(2*(nphi-k+koffset-kghost)+1)*ntheta+(j+joffset-jghost)];
          }
        }
      }
      if (koffset+nphi_proc>=nphi) {
        for(int k = nphi_proc+kghost; k < nphi_proc+2*kghost; k++) {
          int jbeg = (joffset==0) ? jghost : 0;
          int jend = (joffset+ntheta_proc>=ntheta) ? ntheta_proc+jghost : ntheta_proc+2*jghost;
          for(int j = jbeg; j < jend; j++) {
            this->Ylm_r(l,m,k,j) = Yr[(2*(nphi-k+koffset-kghost)+1)*ntheta+(j+joffset-jghost)];
            this->Slm_th(l,m,k,j) = Stheta[(2*(k-nphi+koffset-kghost)+1)*ntheta+(j+joffset-jghost)];
            this->Slm_phi(l,m,k,j) = Sphi[(2*(k-nphi+koffset-kghost)+1)*ntheta+(j+joffset-jghost)];
            this->Tlm_th(l,m,k,j) = Ttheta[(2*(k-nphi+koffset-kghost)+1)*ntheta+(j+joffset-jghost)];
            this->Tlm_phi(l,m,k,j) = Tphi[(2*(k-nphi+koffset-kghost)+1)*ntheta+(j+joffset-jghost)];
          }
        }
      }

      if (joffset==0) {
        real parity = pow(-1,m); // parity of Ylm, 1 when even, -1 if odd, Slm and Tlm have the opposite parity of Ylm (derivative)
        for(int j = 0; j < jghost; j++) {
          int jref = jghost-j+1;
          int kbeg = (koffset==0) ? kghost : 0;
          int kend = (koffset+nphi_proc>=nphi) ? nphi_proc+kghost : nphi_proc+2*kghost;
          for(int k = kbeg; k < kend; k++) {
            this->Ylm_r(l,m,k,j) = parity*Yr[(2*(k+koffset-kghost)+1)*ntheta+(jref+joffset-jghost)];
            this->Slm_th(l,m,k,j) = -parity*Stheta[(2*(k+koffset-kghost)+1)*ntheta+(jref+joffset-jghost)];
            this->Slm_phi(l,m,k,j) = -parity*Sphi[(2*(k+koffset-kghost)+1)*ntheta+(jref+joffset-jghost)];
            this->Tlm_th(l,m,k,j) = -parity*Ttheta[(2*(k+koffset-kghost)+1)*ntheta+(jref+joffset-jghost)];
            this->Tlm_phi(l,m,k,j) = -parity*Tphi[(2*(k+koffset-kghost)+1)*ntheta+(jref+joffset-jghost)];
          }
        }
      }
      if (joffset+ntheta_proc>=ntheta) {
        real parity = pow(-1,m); // parity of Ylm, 1 when even, -1 if odd, Slm and Tlm have the opposite parity of Ylm (derivative)
        for(int j = ntheta_proc+jghost; j < ntheta_proc+2*jghost; j++) {
          int jref = 2*(ntheta_proc+jghost)-(j+1);
          int kbeg = (koffset==0) ? kghost : 0;
          int kend = (koffset+nphi_proc>=nphi) ? nphi_proc+kghost : nphi_proc+2*kghost;
          for(int k = kbeg; k < kend; k++) {
            this->Ylm_r(l,m,k,j) = parity*Yr[(2*(k+koffset-kghost)+1)*ntheta+(jref+joffset-jghost)];
            this->Slm_th(l,m,k,j) = -parity*Stheta[(2*(k+koffset-kghost)+1)*ntheta+(jref+joffset-jghost)];
            this->Slm_phi(l,m,k,j) = -parity*Sphi[(2*(k+koffset-kghost)+1)*ntheta+(jref+joffset-jghost)];
            this->Tlm_th(l,m,k,j) = -parity*Ttheta[(2*(k+koffset-kghost)+1)*ntheta+(jref+joffset-jghost)];
            this->Tlm_phi(l,m,k,j) = -parity*Tphi[(2*(k+koffset-kghost)+1)*ntheta+(jref+joffset-jghost)];
          }
        }
      }
    }
  }
  shtns_free(Yr);
  shtns_free(Tphi);
  shtns_free(Ttheta);
  shtns_free(Sphi);
  shtns_free(Stheta);
  shtns_free(Ylm);
  shtns_free(Tlm);
  shtns_free(Slm);
  shtns_destroy(shtns);
}

void Vsh::GenerateInterfaceGhostVsh() {
  shtns_cfg shtns;                // handle to a sht transform configuration
  int NLM;
  std::complex<real> *Slm, *Tlm;      // spherical harmonics coefficients (l,m space): complex numbers.
  real *Stheta, *Sphi;                // real space : theta,phi
  real *Ttheta, *Tphi;                // real space : theta,phi

  const int mres = 1;             // periodicity in phi (1 for full-sphere, 2 for half the sphere, 3 for 1/3, etc...)

  shtns_verbose(1);                       // displays informations during initialization.
  shtns_use_threads(0);           // enable multi-threaded transforms (if supported).
//        shtns = shtns_init( sht_gauss, lmax, mmax, mres, ntheta, nphi );
//  shtns = shtns_init( sht_reg_fast, lmax, mmax, mres, ntheta_shtns, nphi_shtns );
  shtns = shtns_init( sht_reg_poles, lmax, mmax, mres, ntheta_shtns, nphi_shtns );

  NLM = shtns->nlm;
  Tphi = (real *) shtns_malloc( NSPAT_ALLOC(shtns) * sizeof(real));
  Ttheta = (real *) shtns_malloc( NSPAT_ALLOC(shtns) * sizeof(real));
  Sphi = (real *) shtns_malloc( NSPAT_ALLOC(shtns) * sizeof(real));
  Stheta = (real *) shtns_malloc( NSPAT_ALLOC(shtns) * sizeof(real));
  Slm = (std::complex<real> *) shtns_malloc( NLM * sizeof(std::complex<real>));
  Tlm = (std::complex<real> *) shtns_malloc( NLM * sizeof(std::complex<real>));
  for(int l = 0; l < this->lmax; l++) {
    int sqrt_ll = (l != 0) ? sqrt(l*(l+1)) : 1;
    for(int m = 0; m < this->mmax; m++) {
      LM_LOOP(shtns, Slm[lm]=0.0;  Tlm[lm] = 0.0; )
      Slm[LM(shtns,l,m)] = 1.0/sqrt_ll;
      Tlm[LM(shtns,l,m)] = 1.0/sqrt_ll;
      SHsph_to_spat(shtns,Slm,Stheta,Sphi);
      SHtor_to_spat(shtns,Tlm,Ttheta,Tphi);

      if (koffset==0) {
        for(int k = 0; k < kghost; k++) {
          int jbeg = (joffset==0) ? jghost : 0;
          int jend = (joffset+ntheta_proc>=ntheta) ? ntheta_proc+jghost : ntheta_proc+2*jghost;
          for(int j = jbeg; j < jend; j++) {
            this->Slm_ths(l,m,k,j) = Stheta[(2*(nphi-k+koffset-kghost)+1)*ntheta_shtns+2*(j+joffset-jghost)];
            this->Slm_phis(l,m,k,j) = Sphi[(2*(nphi-k+koffset-kghost))*ntheta_shtns+2*(j+joffset-jghost)+1];
            this->Tlm_ths(l,m,k,j) = Ttheta[(2*(nphi-k+koffset-kghost)+1)*ntheta_shtns+2*(j+joffset-jghost)];
            this->Tlm_phis(l,m,k,j) = Tphi[(2*(nphi-k+koffset-kghost))*ntheta_shtns+2*(j+joffset-jghost)+1];
          }
        }
      }
      if (koffset+nphi_proc>=nphi) {
        for(int k = nphi_proc+kghost; k < nphi_proc+2*kghost; k++) {
          int jbeg = (joffset==0) ? jghost : 0;
          int jend = (joffset+ntheta_proc>=ntheta) ? ntheta_proc+jghost : ntheta_proc+2*jghost;
          for(int j = jbeg; j < jend; j++) {
            this->Slm_ths(l,m,k,j) = Stheta[(2*(k-nphi+koffset-kghost)+1)*ntheta_shtns+2*(j+joffset-jghost)];
            this->Slm_phis(l,m,k,j) = Sphi[(2*(k-nphi+koffset-kghost))*ntheta_shtns+2*(j+joffset-jghost)+1];
            this->Tlm_ths(l,m,k,j) = Ttheta[(2*(k-nphi+koffset-kghost)+1)*ntheta_shtns+2*(j+joffset-jghost)];
            this->Tlm_phis(l,m,k,j) = Tphi[(2*(k-nphi+koffset-kghost))*ntheta_shtns+2*(j+joffset-jghost)+1];
          }
        }
      }

      if (joffset==0) {
        real parity = pow(-1,m); // parity of Ylm, 1 when even, -1 if odd, Slm and Tlm have the opposite parity of Ylm (derivative)
        for(int j = 0; j < jghost; j++) { 
          int jref = jghost-j+1;
          int jrefth = jghost-j+1+1; // empirical, less different cells with T11 it seems
          int kbeg = (koffset==0) ? kghost : 0;
          int kend = (koffset+nphi_proc>=nphi) ? nphi_proc+kghost : nphi_proc+2*kghost;
          for(int k = kbeg; k < kend; k++) {
            this->Slm_ths(l,m,k,j) = -parity*Stheta[(2*(k+koffset-kghost)+1)*ntheta_shtns+2*(jrefth+joffset-jghost)];
            this->Slm_phis(l,m,k,j) = -parity*Sphi[(2*(k+koffset-kghost))*ntheta_shtns+2*(jref+joffset-jghost)+1];
            this->Tlm_ths(l,m,k,j) = -parity*Ttheta[(2*(k+koffset-kghost)+1)*ntheta_shtns+2*(jrefth+joffset-jghost)];
            this->Tlm_phis(l,m,k,j) = -parity*Tphi[(2*(k+koffset-kghost))*ntheta_shtns+2*(jref+joffset-jghost)+1];
          }
        }
      }
      if (joffset+ntheta_proc>=ntheta) {
        for(int j = ntheta_proc+jghost; j < ntheta_proc+2*jghost; j++) {
          real parity = pow(-1,m); // parity of Ylm, 1 when even, -1 if odd, Slm and Tlm have the opposite parity of Ylm (derivative)
          int jref = 2*(ntheta_proc+jghost)-(j+1);
          int jrefth = 2*(ntheta_proc+jghost)-(j+1); // Very very empirical again
          int kbeg = (koffset==0) ? kghost : 0;
          int kend = (koffset+nphi_proc>=nphi) ? nphi_proc+kghost : nphi_proc+2*kghost;
          for(int k = kbeg; k < kend; k++) {
            this->Slm_ths(l,m,k,j) = -parity*Stheta[(2*(k+koffset-kghost)+1)*ntheta_shtns+2*(jrefth+joffset-jghost)];
            this->Slm_phis(l,m,k,j) = -parity*Sphi[(2*(k+koffset-kghost))*ntheta_shtns+2*(jref+joffset-jghost)+1];
            this->Tlm_ths(l,m,k,j) = -parity*Ttheta[(2*(k+koffset-kghost)+1)*ntheta_shtns+2*(jrefth+joffset-jghost)];
            this->Tlm_phis(l,m,k,j) = -parity*Tphi[(2*(k+koffset-kghost))*ntheta_shtns+2*(jref+joffset-jghost)+1];
          }
        }
      }
    }
  }
  shtns_free(Tphi);
  shtns_free(Ttheta);
  shtns_free(Sphi);
  shtns_free(Stheta);
  shtns_free(Tlm);
  shtns_free(Slm);
  shtns_destroy(shtns);
}

void Vsh::write_vect(std::string fn, double *vec, int N)
{
        FILE *fp;
        int i;

        fp = fopen(fn.c_str(),"w");
        for (i=0;i<N;i++) {
                fprintf(fp,"%.6g ",vec[i]);
        }
        fclose(fp);
}

void Vsh::write_mx(std::string fn, double *mx, int N1, int N2)
{
        FILE *fp;
        int i,j;

//        fp = fopen(fn,"w");
        fp = fopen(fn.c_str(),"w");
        for (i=0;i<N1;i++) {
                for(j=0;j<N2;j++) {
                        fprintf(fp,"%.6g ",mx[i*N2+j]);
                }
                fprintf(fp,"\n");
        }
        fclose(fp);
}
