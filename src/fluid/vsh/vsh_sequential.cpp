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

VSH::VSH() {
// Default constructor
} 

void VSH::InitVSH(int nphi, int ntheta, int koffset, int joffset, int kghost, int jghost, int lmax, int mmax) {
  this->nphi = nphi;
  this->nphi_shtns = 2*nphi; //so to have the points at the right coordinates, because the phi-grid is equally spaced between 0 (included) and 2*M_PI (excluded)
  this->ntheta = ntheta;
  this->ntheta_shtns = 2*ntheta + 1;
  this->koffset = koffset;
  this->joffset = joffset;
  this->kghost = kghost;
  this->jghost = jghost;
  this->lmax = lmax;
  this->mmax = mmax;

  this->Ylm_r = IdefixArray4D<real> ("Ylm_r", lmax, mmax, nphi+2*kghost, ntheta+2*jghost);
  this->Slm_th = IdefixArray4D<real> ("Slm_th", lmax, mmax, nphi+2*kghost, ntheta+2*jghost);
  this->Slm_phi = IdefixArray4D<real> ("Slm_phi", lmax, mmax, nphi+2*kghost, ntheta+2*jghost);
  this->Tlm_th = IdefixArray4D<real> ("Tlm_th", lmax, mmax, nphi+2*kghost, ntheta+2*jghost);
  this->Tlm_phi = IdefixArray4D<real> ("Tlm_phi", lmax, mmax, nphi+2*kghost, ntheta+2*jghost);

  this->Slm_ths = IdefixArray4D<real> ("Slm_ths", lmax, mmax, nphi+2*kghost, ntheta+2*jghost);
  this->Slm_phis = IdefixArray4D<real> ("Slm_phis", lmax, mmax, nphi+2*kghost, ntheta+2*jghost);
  this->Tlm_ths = IdefixArray4D<real> ("Tlm_ths", lmax, mmax, nphi+2*kghost, ntheta+2*jghost);
  this->Tlm_phis = IdefixArray4D<real> ("Tlm_phis", lmax, mmax, nphi+2*kghost, ntheta+2*jghost);
}

void VSH::GenerateCellVSH(int write) {
  shtns_cfg shtns;                // handle to a sht transform configuration
  int NLM;
  std::complex<real> *Ylm, *Slm, *Tlm;      // spherical harmonics coefficients (l,m space): complex numbers.
  real *Yr;                // real space : r
  real *Stheta, *Sphi;                // real space : theta,phi
  real *Ttheta, *Tphi;                // real space : theta,phi
  real t;

  const int mres = 1;             // periodicity in phi (1 for full-sphere, 2 for half the sphere, 3 for 1/3, etc...)

  shtns_verbose(1);                       // displays informations during initialization.
  shtns_use_threads(0);           // enable multi-threaded transforms (if supported).
//        shtns = shtns_init( sht_gauss, lmax, mmax, mres, ntheta, nphi );
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
    for(int m = 0; m < this->mmax; m++) {
      LM_LOOP(shtns,  Ylm[lm]=0.0; Slm[lm]=0.0;  Tlm[lm] = 0.0; )
      Ylm[LM(shtns,l,m)] = 1.0;
      Slm[LM(shtns,l,m)] = 1.0;
      Tlm[LM(shtns,l,m)] = 1.0;
      SH_to_spat(shtns,Ylm,Yr);
      SHsph_to_spat(shtns,Slm,Stheta,Sphi);
      SHtor_to_spat(shtns,Tlm,Ttheta,Tphi);
      if (write==1) {
        std::string path;
        path = "datavsh/Yr"+std::to_string(l)+std::to_string(m);
        std::cout << path << std::endl;
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
      for(int k = kghost; k < nphi+kghost; k++) {
        for(int j = jghost; j < ntheta+jghost; j++) {
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

void VSH::GenerateInterfaceVSH(int write) { //defined at the right interface between cells in the direction of the corresponding component, like the magnetic field BX1s,2s,3s
  shtns_cfg shtns;                // handle to a sht transform configuration
  int NLM;
  std::complex<real> *Slm, *Tlm;      // spherical harmonics coefficients (l,m space): complex numbers.
  real *Stheta, *Sphi;                // real space : theta,phi
  real *Ttheta, *Tphi;                // real space : theta,phi
  real t;

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
    for(int m = 0; m < this->mmax; m++) {
      LM_LOOP(shtns, Slm[lm]=0.0;  Tlm[lm] = 0.0; )
      Slm[LM(shtns,l,m)] = 1.0;
      Tlm[LM(shtns,l,m)] = 1.0;
      SHsph_to_spat(shtns,Slm,Stheta,Sphi);
      SHtor_to_spat(shtns,Tlm,Ttheta,Tphi);
      if (write==1) {
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
      for(int k = kghost; k < nphi+kghost; k++) {
        for(int j = jghost; j < ntheta+jghost; j++) {
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

void VSH::GenerateCellGhostVSH() {
  for(int l = 0; l < this->lmax; l++) {
    for(int m = 0; m < this->mmax; m++) {
      for(int k = 0; k < kghost; k++) {
        for(int j = jghost; j < ntheta+jghost; j++) {
          this->Ylm_r(l,m,k,j) = this->Ylm_r(l,m,nphi-k,j);
          this->Slm_th(l,m,k,j) = this->Slm_th(l,m,nphi-k,j);
          this->Slm_phi(l,m,k,j) = this->Slm_phi(l,m,nphi-k,j);
          this->Tlm_th(l,m,k,j) = this->Tlm_th(l,m,nphi-k,j);
          this->Tlm_phi(l,m,k,j) = this->Tlm_phi(l,m,nphi-k,j);
        }
      }
      for(int k = nphi+kghost; k < nphi+2*kghost; k++) {
        for(int j = jghost; j < ntheta+jghost; j++) {
          this->Ylm_r(l,m,k,j) = this->Ylm_r(l,m,k-nphi,j);
          this->Slm_th(l,m,k,j) = this->Slm_th(l,m,k-nphi,j);
          this->Slm_phi(l,m,k,j) = this->Slm_phi(l,m,k-nphi,j);
          this->Tlm_th(l,m,k,j) = this->Tlm_th(l,m,k-nphi,j);
          this->Tlm_phi(l,m,k,j) = this->Tlm_phi(l,m,k-nphi,j);
        }
      }

      real parity = pow(-1,m); // parity of Ylm, 1 when even, -1 if odd, Slm and Tlm have the opposite parity of Ylm (derivative)
      for(int j = 0; j < jghost; j++) {
        int jref = jghost-j+1;
        for(int k = kghost; k < nphi+kghost; k++) {
          this->Ylm_r(l,m,k,j) = parity*this->Ylm_r(l,m,k,jref);
          this->Slm_th(l,m,k,j) = -parity*this->Slm_th(l,m,k,jref);
          this->Slm_phi(l,m,k,j) = -parity*this->Slm_phi(l,m,k,jref);
          this->Tlm_th(l,m,k,j) = -parity*this->Tlm_th(l,m,k,jref);
          this->Tlm_phi(l,m,k,j) = -parity*this->Tlm_phi(l,m,k,jref);
        }
      }
      for(int j = ntheta+jghost; j < ntheta+2*jghost; j++) {
        int jref = 2*(ntheta+jghost)-(j+1);
        for(int k = kghost; k < nphi+kghost; k++) {
          this->Ylm_r(l,m,k,j) = parity*this->Ylm_r(l,m,k,jref);
          this->Slm_th(l,m,k,j) = -parity*this->Slm_th(l,m,k,jref);
          this->Slm_phi(l,m,k,j) = -parity*this->Slm_phi(l,m,k,jref);
          this->Tlm_th(l,m,k,j) = -parity*this->Tlm_th(l,m,k,jref);
          this->Tlm_phi(l,m,k,j) = -parity*this->Tlm_phi(l,m,k,jref);
        }
      }
    }
  }
}

void VSH::GenerateInterfaceGhostVSH() {
  for(int l = 0; l < this->lmax; l++) {
    for(int m = 0; m < this->mmax; m++) {
      for(int k = 0; k < kghost; k++) {
        for(int j = jghost; j < ntheta+jghost; j++) {
          this->Slm_ths(l,m,k,j) = this->Slm_ths(l,m,nphi-k,j);
          this->Slm_phis(l,m,k,j) = this->Slm_phis(l,m,nphi-k,j);
          this->Tlm_ths(l,m,k,j) = this->Tlm_ths(l,m,nphi-k,j);
          this->Tlm_phis(l,m,k,j) = this->Tlm_phis(l,m,nphi-k,j);
        }
      }
      for(int k = nphi+kghost; k < nphi+2*kghost; k++) {
        for(int j = jghost; j < ntheta+jghost; j++) {
          this->Slm_ths(l,m,k,j) = this->Slm_ths(l,m,k-nphi,j);
          this->Slm_phis(l,m,k,j) = this->Slm_phis(l,m,k-nphi,j);
          this->Tlm_ths(l,m,k,j) = this->Tlm_ths(l,m,k-nphi,j);
          this->Tlm_phis(l,m,k,j) = this->Tlm_phis(l,m,k-nphi,j);
        }
      }

      real parity = pow(-1,m); // parity of Ylm, 1 when even, -1 if odd, Slm and Tlm have the opposite parity of Ylm (derivative)
      for(int j = 0; j < jghost; j++) { 
        int jref = jghost-j+1;
        int jrefth = jghost-j+1+1;
//        int jrefth = jghost-j+1;
        for(int k = kghost; k < nphi+kghost; k++) {
          this->Slm_ths(l,m,k,j) = -parity*this->Slm_ths(l,m,k,jrefth);
          this->Slm_phis(l,m,k,j) = -parity*this->Slm_phis(l,m,k,jref);
          this->Tlm_ths(l,m,k,j) = -parity*this->Tlm_ths(l,m,k,jrefth);
          this->Tlm_phis(l,m,k,j) = -parity*this->Tlm_phis(l,m,k,jref);
        }
      }
      for(int j = ntheta+jghost; j < ntheta+2*jghost; j++) {
        int jref = 2*(ntheta+jghost)-(j+1);
        int jrefth = 2*(ntheta+jghost)-(j+1)+1;
//        int jrefth = 2*(ntheta+jghost)-(j+1);
        for(int k = kghost; k < nphi+kghost; k++) {
          this->Slm_ths(l,m,k,j) = -parity*this->Slm_ths(l,m,k,jrefth);
          this->Slm_phis(l,m,k,j) = -parity*this->Slm_phis(l,m,k,jref);
          this->Tlm_ths(l,m,k,j) = -parity*this->Tlm_ths(l,m,k,jrefth);
          this->Tlm_phis(l,m,k,j) = -parity*this->Tlm_phis(l,m,k,jref);
        }
      }
    }
  }
}

void VSH::write_vect(std::string fn, double *vec, int N)
{
        FILE *fp;
        int i;

        fp = fopen(fn.c_str(),"w");
        for (i=0;i<N;i++) {
                fprintf(fp,"%.6g ",vec[i]);
        }
        fclose(fp);
}

void VSH::write_mx(std::string fn, double *mx, int N1, int N2)
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
