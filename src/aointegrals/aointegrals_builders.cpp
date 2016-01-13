/* 
 *  The Chronus Quantum (ChronusQ) software package is high-performace 
 *  computational chemistry software with a strong emphasis on explicitly 
 *  time-dependent and post-SCF quantum mechanical methods.
 *  
 *  Copyright (C) 2014-2015 Li Research Group (University of Washington)
 *  
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *  
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *  
 *  You should have received a copy of the GNU General Public License along
 *  with this program; if not, write to the Free Software Foundation, Inc.,
 *  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *  
 *  Contact the Developers:
 *    E-Mail: xsli@uw.edu
 *  
 */
#include <aointegrals.h>
using ChronusQ::AOIntegrals;

void AOIntegrals::OneEDriver(OneBodyEngine::integral_type iType) {
  std::vector<RealMap> mat;
  int NB = this->nTCS_*this->nBasis_;
  int NBSq = NB*NB;
  if(iType == OneBodyEngine::overlap){
    mat.push_back(RealMap(this->overlap_->data(),NB,NB));
  } else if(iType == OneBodyEngine::kinetic) {
    mat.push_back(RealMap(this->kinetic_->data(),NB,NB));
  } else if(iType == OneBodyEngine::nuclear) {
    mat.push_back(RealMap(this->potential_->data(),NB,NB));
  } else if(iType == OneBodyEngine::emultipole1) {
    mat.push_back(RealMap(this->overlap_->data(),NB,NB));
    for(auto i = 0, IOff=0; i < 3; i++,IOff+=NBSq)
      mat.push_back(RealMap(&this->elecDipole_->storage()[IOff],NB,NB));
  } else if(iType == OneBodyEngine::emultipole2) {
    mat.push_back(RealMap(this->overlap_->data(),NB,NB));
    for(auto i = 0, IOff=0; i < 3; i++,IOff+=NBSq)
      mat.push_back(RealMap(&this->elecDipole_->storage()[IOff],NB,NB));
    for(auto i = 0, IOff=0; i < 6; i++,IOff+=NBSq)
      mat.push_back(RealMap(&this->elecQuadpole_->storage()[IOff],NB,NB));
  } else if(iType == OneBodyEngine::emultipole3) {
    mat.push_back(RealMap(this->overlap_->data(),NB,NB));
    for(auto i = 0, IOff=0; i < 3; i++,IOff+=NBSq)
      mat.push_back(RealMap(&this->elecDipole_->storage()[IOff],NB,NB));
    for(auto i = 0, IOff=0; i < 6; i++,IOff+=NBSq)
      mat.push_back(RealMap(&this->elecQuadpole_->storage()[IOff],NB,NB));
    for(auto i = 0, IOff=0; i < 10; i++,IOff+=NBSq)
      mat.push_back(RealMap(&this->elecOctpole_->storage()[IOff],NB,NB));
  } else {
    cout << "OneBodyEngine type not recognized" << endl;
    exit(EXIT_FAILURE);
  }
 
#ifdef _OPENMP
  int nthreads = omp_get_max_threads();
#else
  int nthreads = 1;
#endif
  // Define integral Engine
  std::vector<OneBodyEngine> engines(nthreads);
  engines[0] = OneBodyEngine(iType,this->basisSet_->maxPrim(),this->basisSet_->maxL(),0);

  // If engine is V, define nuclear charges
  if(iType == OneBodyEngine::nuclear){
    std::vector<std::pair<double,std::array<double,3>>> q;
    for(int i = 0; i < this->molecule_->nAtoms(); i++) {
      q.push_back(
        {
          static_cast<double>((*this->molecularConstants_).atomZ[i]), 
          {
            {
	      (*this->molecularConstants_).cart[0][i],
	      (*this->molecularConstants_).cart[1][i],
	      (*this->molecularConstants_).cart[2][i]
	    }
	  }
	}
      );
    }
    engines[0].set_q(q);
  }
  for(size_t i = 1; i < nthreads; i++) engines[i] = engines[0];

  if(!this->basisSet_->haveMapSh2Bf) this->basisSet_->makeMapSh2Bf(this->nTCS_); 
#ifdef _OPENMP
  #pragma omp parallel
#endif
  {
#ifdef _OPENMP
    int thread_id = omp_get_thread_num();
#else
    int thread_id = 0;
#endif
    for(auto s1=0l, s12=0l; s1 < this->basisSet_->nShell(); s1++){
      int bf1_s = this->basisSet_->mapSh2Bf(s1);
      int n1  = this->basisSet_->shells(s1).size();
      for(int s2=0; s2 <= s1; s2++, s12++){
        if(s12 % nthreads != thread_id) continue;
        int bf2_s = this->basisSet_->mapSh2Bf(s2);
        int n2  = this->basisSet_->shells(s2).size();
  
        const double* buff = engines[thread_id].compute(
          this->basisSet_->shells(s1),
          this->basisSet_->shells(s2)
        );

        int IOff = 0;
        for(auto nMat = 0; nMat < mat.size(); nMat++) {
//        ConstRealMap bufMat(&buff[IOff],n1,n2); // Read only map
          Eigen::Map<const Eigen::Matrix<double,Dynamic,Dynamic,Eigen::RowMajor>>
            bufMat(&buff[IOff],n1,n2);
          for(auto i = 0, bf1 = bf1_s; i < n1; i++, bf1 += this->nTCS_)            
          for(auto j = 0, bf2 = bf2_s; j < n2; j++, bf2 += this->nTCS_){            
            mat[nMat](bf1,bf2) = bufMat(i,j);
            if(this->nTCS_ == 2) mat[nMat](bf1+1,bf2+1) = bufMat(i,j);
          }
          IOff += n1*n2;
        }
      }
    }
  } // end openmp parallel
  for(auto nMat = 0; nMat < mat.size(); nMat++) 
    mat[nMat] = mat[nMat].selfadjointView<Lower>();
}

void AOIntegrals::computeAOOneE(){
  // Collect Relevant data into a struct (odd, but convienient) 
  this->iniMolecularConstants();
/*
  BasisSet newBasis;
  newBasis.communicate(*this->fileio_);
  this->basisSet_->genUCvomLocal(&newBasis);
  newBasis.makeMaps(this->nTCS_,this->molecule_);
  if(this->isPrimary) {
    cout << "Old Basis" << endl;
    for(auto i = 0 ; i < this->basisSet_->nShell(); i++)
      cout << this->basisSet_->shells(i) << endl;
    cout << "New Basis" << endl;
    for(auto i = 0 ; i < newBasis.nShell(); i++)
      cout << newBasis.shells(i) << endl;
  }
*/

  // Start timer for one-electron integral evaluation
  auto oneEStart = std::chrono::high_resolution_clock::now();

  // Compute and time overlap integrals
  auto OStart = std::chrono::high_resolution_clock::now();
  if(this->maxMultipole_ ==3) OneEDriver(OneBodyEngine::emultipole3);
  else if(this->maxMultipole_ == 2) OneEDriver(OneBodyEngine::emultipole2);
  else if(this->maxMultipole_ == 1) OneEDriver(OneBodyEngine::emultipole1);
  else OneEDriver(OneBodyEngine::overlap);
  auto OEnd = std::chrono::high_resolution_clock::now();

  // Compute and time kinetic integrals
  auto TStart = std::chrono::high_resolution_clock::now();
  OneEDriver(OneBodyEngine::kinetic);
  auto TEnd = std::chrono::high_resolution_clock::now();

  // Compute and time nuclear attraction integrals (negative sign is factored in)
  auto VStart = std::chrono::high_resolution_clock::now();
  OneEDriver(OneBodyEngine::nuclear);
  auto VEnd = std::chrono::high_resolution_clock::now();

// add DKH correction to kinetic energy
//  if (this->isPrimary) this->DKH0();

// Build Core Hamiltonian
  (*this->oneE_) = (*this->kinetic_) + (*this->potential_);

  // Get end time of one-electron integral evaluation
  auto oneEEnd = std::chrono::high_resolution_clock::now();
  if(this->controls_->printLevel >= 2) this->printOneE();

  // Compute time differenes
  this->OneED = oneEEnd - oneEStart;
  this->SED = OEnd - OStart;
  this->TED = TEnd - TStart;
  this->VED = VEnd - VStart;
  this->haveAOOneE = true;
  if(this->isPrimary) this->writeOneE();
}

using libint2::TwoBodyEngine;
void AOIntegrals::computeSchwartz(){
  if(getRank() == 0) {
    RealMatrix *ShBlk; 
    this->schwartz_->setZero();
 
    // Define Integral Engine
    TwoBodyEngine<libint2::Coulomb> engine = 
      TwoBodyEngine<libint2::Coulomb>(this->basisSet_->maxPrim(),
                                      this->basisSet_->maxL(),0);
    engine.set_precision(0.); // Don't screen primitives during schwartz
 
    auto start =  std::chrono::high_resolution_clock::now();
    for(int s1=0; s1 < this->basisSet_->nShell(); s1++){
      int n1  = this->basisSet_->shells(s1).size();
      for(int s2=0; s2 <= s1; s2++){
        int n2  = this->basisSet_->shells(s2).size();
 
        const auto* buff = engine.compute(
          this->basisSet_->shells(s1),
          this->basisSet_->shells(s2),
          this->basisSet_->shells(s1),
          this->basisSet_->shells(s2)
        );
 
        
        ShBlk = new RealMatrix(n1,n2);
        ShBlk->setZero();
 
        for(auto i = 0, ij = 0; i < n1; i++)
        for(auto j = 0; j < n2; j++, ij++){
          (*ShBlk)(i,j) = buff[ij*n1*n2 + ij];
        }
 
        (*this->schwartz_)(s1,s2) = std::sqrt(ShBlk->lpNorm<Infinity>());
        
        delete ShBlk;
      }
    }
    auto finish =  std::chrono::high_resolution_clock::now();
    this->SchwartzD = finish - start;
    (*this->schwartz_) = this->schwartz_->selfadjointView<Lower>();
  }
#ifdef CQ_ENABLE_MPI
  MPI_Bcast(this->schwartz_->data(),this->schwartz_->size(),MPI_DOUBLE,0,
    MPI_COMM_WORLD);
#endif

  this->haveSchwartz = true;
}
void AOIntegrals::computeAOTwoE(){
  if(!this->haveSchwartz) this->computeSchwartz();
  if(getRank() != 0) return;

#ifdef _OPENMP
  int nthreads = omp_get_max_threads();
#else
  int nthreads = 1;
#endif
  std::vector<coulombEngine> engines(nthreads);
  engines[0] = coulombEngine(this->basisSet_->maxPrim(),this->basisSet_->maxL(),0);
  engines[0].set_precision(std::numeric_limits<double>::epsilon());

  for(int i=1; i<nthreads; i++) engines[i] = engines[0];
  if(!this->basisSet_->haveMapSh2Bf) this->basisSet_->makeMapSh2Bf(this->nTCS_); 

  this->aoERI_->fill(0.0);

#ifdef _OPENMP
  #pragma omp parallel
#endif
  {
#ifdef _OPENMP
    int thread_id = omp_get_thread_num();
#else
    int thread_id = 0;
#endif
  for(int s1 = 0, s1234=0; s1 < this->basisSet_->nShell(); s1++) {
    int bf1_s = this->basisSet_->mapSh2Bf(s1);
    int n1    = this->basisSet_->shells(s1).size();
    for(int s2 = 0; s2 <= s1; s2++) {
      int bf2_s = this->basisSet_->mapSh2Bf(s2);
      int n2    = this->basisSet_->shells(s2).size();
      for(int s3 = 0; s3 <= s1; s3++) {
        int bf3_s = this->basisSet_->mapSh2Bf(s3);
        int n3    = this->basisSet_->shells(s3).size();
        int s4_max = (s1 == s3) ? s2 : s3;
        for(int s4 = 0; s4 <= s4_max; s4++, s1234++) {

          if(s1234 % nthreads != thread_id) continue;

          int bf4_s = this->basisSet_->mapSh2Bf(s4);
          int n4    = this->basisSet_->shells(s4).size();
    
          // Schwartz and Density screening
          if((*this->schwartz_)(s1,s2) * (*this->schwartz_)(s3,s4)
              < this->controls_->thresholdSchawrtz ) continue;
 
          const double* buff = engines[thread_id].compute(
            this->basisSet_->shells(s1),
            this->basisSet_->shells(s2),
            this->basisSet_->shells(s3),
            this->basisSet_->shells(s4));


          std::vector<std::array<int,4>> lower;
          std::vector<std::array<int,4>> upper;
          for(int i = 0, bf1 = bf1_s, ijkl = 0 ; i < n1; ++i, bf1 += this->nTCS_) 
          for(int j = 0, bf2 = bf2_s           ; j < n2; ++j, bf2 += this->nTCS_) 
          for(int k = 0, bf3 = bf3_s           ; k < n3; ++k, bf3 += this->nTCS_) 
          for(int l = 0, bf4 = bf4_s           ; l < n4; ++l, bf4 += this->nTCS_, ++ijkl) {
            (*this->aoERI_)(bf1,bf2,bf3,bf4) = buff[ijkl];
            (*this->aoERI_)(bf1,bf2,bf4,bf3) = buff[ijkl];
            (*this->aoERI_)(bf2,bf1,bf3,bf4) = buff[ijkl];
            (*this->aoERI_)(bf2,bf1,bf4,bf3) = buff[ijkl];
            (*this->aoERI_)(bf3,bf4,bf1,bf2) = buff[ijkl];
            (*this->aoERI_)(bf4,bf3,bf1,bf2) = buff[ijkl];
            (*this->aoERI_)(bf3,bf4,bf2,bf1) = buff[ijkl];
            (*this->aoERI_)(bf4,bf3,bf2,bf1) = buff[ijkl];
            if(this->nTCS_ == 2){
              (*this->aoERI_)(bf1+1,bf2+1,bf3+1,bf4+1) = buff[ijkl];
              (*this->aoERI_)(bf1+1,bf2+1,bf4+1,bf3+1) = buff[ijkl];
              (*this->aoERI_)(bf2+1,bf1+1,bf3+1,bf4+1) = buff[ijkl];
              (*this->aoERI_)(bf2+1,bf1+1,bf4+1,bf3+1) = buff[ijkl];
              (*this->aoERI_)(bf3+1,bf4+1,bf1+1,bf2+1) = buff[ijkl];
              (*this->aoERI_)(bf4+1,bf3+1,bf1+1,bf2+1) = buff[ijkl];
              (*this->aoERI_)(bf3+1,bf4+1,bf2+1,bf1+1) = buff[ijkl];
              (*this->aoERI_)(bf4+1,bf3+1,bf2+1,bf1+1) = buff[ijkl];

              (*this->aoERI_)(bf1+1,bf2+1,bf3,bf4)     = buff[ijkl];
              (*this->aoERI_)(bf1+1,bf2+1,bf4,bf3)     = buff[ijkl];
              (*this->aoERI_)(bf2+1,bf1+1,bf3,bf4)     = buff[ijkl];
              (*this->aoERI_)(bf2+1,bf1+1,bf4,bf3)     = buff[ijkl];
              (*this->aoERI_)(bf3+1,bf4+1,bf1,bf2)     = buff[ijkl];
              (*this->aoERI_)(bf4+1,bf3+1,bf1,bf2)     = buff[ijkl];
              (*this->aoERI_)(bf3+1,bf4+1,bf2,bf1)     = buff[ijkl];
              (*this->aoERI_)(bf4+1,bf3+1,bf2,bf1)     = buff[ijkl];

              (*this->aoERI_)(bf1,bf2,bf3+1,bf4+1)     = buff[ijkl];
              (*this->aoERI_)(bf1,bf2,bf4+1,bf3+1)     = buff[ijkl];
              (*this->aoERI_)(bf2,bf1,bf3+1,bf4+1)     = buff[ijkl];
              (*this->aoERI_)(bf2,bf1,bf4+1,bf3+1)     = buff[ijkl];
              (*this->aoERI_)(bf3,bf4,bf1+1,bf2+1)     = buff[ijkl];
              (*this->aoERI_)(bf4,bf3,bf1+1,bf2+1)     = buff[ijkl];
              (*this->aoERI_)(bf3,bf4,bf2+1,bf1+1)     = buff[ijkl];
              (*this->aoERI_)(bf4,bf3,bf2+1,bf1+1)     = buff[ijkl];
            }
	  }

	}
      }
    }
  }
  } // OMP Parallel
  this->haveAOTwoE = true;

//this->fileio_->out << "Two-Electron Integrals (ERIs)" << endl;
//for(auto i = 0; i < this->nBasis_; i++)
//for(auto j = 0; j < this->nBasis_; j++)
//for(auto k = 0; k < this->nBasis_; k++)
//for(auto l = 0; l < this->nBasis_; l++){
//  this->fileio_->out << "(" << i << "," << j << "|" << k << "," << l << ")  ";
//  this->fileio_->out << (*this->aoERI_)(i,j,k,l) << endl;
//};
}



