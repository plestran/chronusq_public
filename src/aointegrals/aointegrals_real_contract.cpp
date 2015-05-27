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
namespace ChronusQ{
#ifdef USE_LIBINT
template<>
void AOIntegrals::twoEContractDirect(bool doRHFFock, const RealMatrix &X, RealMatrix &AX) {
  this->fileio_->out << "Contracting Directly with two-electron integrals" << endl;

  if(!this->haveSchwartz) this->computeSchwartz();
  if(!this->basisSet_->haveMap) this->basisSet_->makeMap(this->molecule_); 
  AX.setZero();
  std::vector<RealMatrix> 
    G(this->controls_->nthreads,RealMatrix::Zero(this->nBasis_,this->nBasis_));

  std::vector<coulombEngine> engines(this->controls_->nthreads);
  engines[0] = coulombEngine(this->basisSet_->maxPrim,this->basisSet_->maxL,0);
  engines[0].set_precision(std::numeric_limits<double>::epsilon());

  for(int i=1; i<this->controls_->nthreads; i++) engines[i] = engines[0];

  auto start = std::chrono::high_resolution_clock::now();
  this->basisSet_->computeShBlkNorm(this->molecule_,&X);
  auto finish = std::chrono::high_resolution_clock::now();
  if(doRHFFock) this->DenShBlkD = finish - start;
  int ijkl = 0;
  start = std::chrono::high_resolution_clock::now();

  auto efficient_twoe = [&] (int thread_id) {
    coulombEngine &engine = engines[thread_id];
    RealMatrix &g = G[thread_id];
    for(int s1 = 0, s1234=0; s1 < this->basisSet_->nShell(); s1++) {
      int bf1_s = this->basisSet_->mapSh2Bf[s1];
      int n1    = this->basisSet_->shells_libint[s1].size();
      for(int s2 = 0; s2 <= s1; s2++) {
        int bf2_s = this->basisSet_->mapSh2Bf[s2];
        int n2    = this->basisSet_->shells_libint[s2].size();
        for(int s3 = 0; s3 <= s1; s3++) {
          int bf3_s = this->basisSet_->mapSh2Bf[s3];
          int n3    = this->basisSet_->shells_libint[s3].size();
          int s4_max = (s1 == s3) ? s2 : s3;
          for(int s4 = 0; s4 <= s4_max; s4++, s1234++) {
            if(s1234 % this->controls_->nthreads != thread_id) continue;
            int bf4_s = this->basisSet_->mapSh2Bf[s4];
            int n4    = this->basisSet_->shells_libint[s4].size();
      
            // Schwartz and Density screening
            if( std::max((*this->basisSet_->shBlkNorm)(s1,s4),
                   std::max((*this->basisSet_->shBlkNorm)(s2,s4),
                      std::max((*this->basisSet_->shBlkNorm)(s3,s4),
                         std::max((*this->basisSet_->shBlkNorm)(s1,s3),
                            std::max((*this->basisSet_->shBlkNorm)(s2,s3),
                                     (*this->basisSet_->shBlkNorm)(s1,s2))
                            )
                         )      
                      )
                   ) * (*this->schwartz_)(s1,s2)
                     * (*this->schwartz_)(s3,s4)
                   < this->controls_->thresholdSchawrtz ) continue;
 
            const double* buff = engine.compute(
              this->basisSet_->shells_libint[s1],
              this->basisSet_->shells_libint[s2],
              this->basisSet_->shells_libint[s3],
              this->basisSet_->shells_libint[s4]);
      
            double s12_deg = (s1 == s2) ? 1.0 : 2.0;
            double s34_deg = (s3 == s4) ? 1.0 : 2.0;
            double s12_34_deg = (s1 == s3) ? (s2 == s4 ? 1.0 : 2.0) : 2.0;
            double s1234_deg = s12_deg * s34_deg * s12_34_deg;
            for(int i = 0, ijkl = 0 ; i < n1; ++i) {
              int bf1 = bf1_s + i;
              for(int j = 0; j < n2; ++j) {
                int bf2 = bf2_s + j;
                for(int k = 0; k < n3; ++k) {
                  int bf3 = bf3_s + k;
                  for(int l = 0; l < n4; ++l, ++ijkl) {
                    int bf4 = bf4_s + l;
                    double v = buff[ijkl]*s1234_deg;

                    // Coulomb
                    g(bf1,bf2) += X(bf3,bf4)*v;
                    g(bf3,bf4) += X(bf1,bf2)*v;
                    g(bf2,bf1) += X(bf4,bf3)*v;
                    g(bf4,bf3) += X(bf2,bf1)*v;
 
                    // Exchange
                    if(doRHFFock) {
                      g(bf1,bf3) -= 0.25*X(bf2,bf4)*v;
                      g(bf4,bf2) -= 0.25*X(bf1,bf3)*v;
                      g(bf1,bf4) -= 0.25*X(bf2,bf3)*v;
                      g(bf3,bf2) -= 0.25*X(bf1,bf4)*v;

                      g(bf3,bf1) -= 0.25*X(bf4,bf2)*v;
                      g(bf2,bf4) -= 0.25*X(bf3,bf1)*v;
                      g(bf4,bf1) -= 0.25*X(bf3,bf2)*v;
                      g(bf2,bf3) -= 0.25*X(bf4,bf1)*v;
                    } else {
                      g(bf1,bf3) -= 0.5*X(bf2,bf4)*v;
                      g(bf4,bf2) -= 0.5*X(bf1,bf3)*v;
                      g(bf1,bf4) -= 0.5*X(bf2,bf3)*v;
                      g(bf3,bf2) -= 0.5*X(bf1,bf4)*v;

                      g(bf3,bf1) -= 0.5*X(bf4,bf2)*v;
                      g(bf2,bf4) -= 0.5*X(bf3,bf1)*v;
                      g(bf4,bf1) -= 0.5*X(bf3,bf2)*v;
                      g(bf2,bf3) -= 0.5*X(bf4,bf1)*v;
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  };

  auto simple_twoe = [&] (int thread_id){
    coulombEngine &engine = engines[thread_id];
    RealMatrix &g = G[thread_id];
    for(int s1 = 0, s1234=0; s1 < this->basisSet_->nShell(); s1++) {
      int bf1_s = this->basisSet_->mapSh2Bf[s1];
      int n1    = this->basisSet_->shells_libint[s1].size();
      for(int s2 = 0; s2 < this->basisSet_->nShell(); s2++) {
        int bf2_s = this->basisSet_->mapSh2Bf[s2];
        int n2    = this->basisSet_->shells_libint[s2].size();
        for(int s3 = 0; s3 < this->basisSet_->nShell(); s3++) {
          int bf3_s = this->basisSet_->mapSh2Bf[s3];
          int n3    = this->basisSet_->shells_libint[s3].size();
          for(int s4 = 0; s4 < this->basisSet_->nShell(); s4++, s1234++) {
            if(s1234 % this->controls_->nthreads != thread_id) continue;
            int bf4_s = this->basisSet_->mapSh2Bf[s4];
            int n4    = this->basisSet_->shells_libint[s4].size();
      
            // Schwartz and Density screening
            if( std::max((*this->basisSet_->shBlkNorm)(s1,s4),
                   std::max((*this->basisSet_->shBlkNorm)(s2,s4),
                      std::max((*this->basisSet_->shBlkNorm)(s3,s4),
                         std::max((*this->basisSet_->shBlkNorm)(s1,s3),
                            std::max((*this->basisSet_->shBlkNorm)(s2,s3),
                                     (*this->basisSet_->shBlkNorm)(s1,s2))
                            )
                         )      
                      )
                   ) * (*this->schwartz_)(s1,s2)
                     * (*this->schwartz_)(s3,s4)
                   < this->controls_->thresholdSchawrtz ) continue;
 
            const double* buff = engine.compute(
              this->basisSet_->shells_libint[s1],
              this->basisSet_->shells_libint[s2],
              this->basisSet_->shells_libint[s3],
              this->basisSet_->shells_libint[s4]);
            for(int i = 0, ijkl = 0 ; i < n1; ++i) {
              int bf1 = bf1_s + i;
              for(int j = 0; j < n2; ++j) {
                int bf2 = bf2_s + j;
                for(int k = 0; k < n3; ++k) {
                  int bf3 = bf3_s + k;
                  for(int l = 0; l < n4; ++l, ++ijkl) {
                    int bf4 = bf4_s + l;
                    double v = buff[ijkl];

                    // Coulomb
                    g(bf1,bf2) += X(bf3,bf4)*v;
                    g(bf3,bf4) += X(bf1,bf2)*v;
                    g(bf2,bf1) += X(bf4,bf3)*v;
                    g(bf4,bf3) += X(bf2,bf1)*v;
 
                    // Exchange
                    if(doRHFFock) {
                      g(bf1,bf3) -= 0.25*X(bf2,bf4)*v;
                      g(bf4,bf2) -= 0.25*X(bf1,bf3)*v;
                      g(bf1,bf4) -= 0.25*X(bf2,bf3)*v;
                      g(bf3,bf2) -= 0.25*X(bf1,bf4)*v;

                      g(bf3,bf1) -= 0.25*X(bf4,bf2)*v;
                      g(bf2,bf4) -= 0.25*X(bf3,bf1)*v;
                      g(bf4,bf1) -= 0.25*X(bf3,bf2)*v;
                      g(bf2,bf3) -= 0.25*X(bf4,bf1)*v;
                    } else {
                      g(bf1,bf3) -= 0.5*X(bf2,bf4)*v;
                      g(bf4,bf2) -= 0.5*X(bf1,bf3)*v;
                      g(bf1,bf4) -= 0.5*X(bf2,bf3)*v;
                      g(bf3,bf2) -= 0.5*X(bf1,bf4)*v;

                      g(bf3,bf1) -= 0.5*X(bf4,bf2)*v;
                      g(bf2,bf4) -= 0.5*X(bf3,bf1)*v;
                      g(bf4,bf1) -= 0.5*X(bf3,bf2)*v;
                      g(bf2,bf3) -= 0.5*X(bf4,bf1)*v;
                    }
                    
                  }
                }
              }
            }

          }
        }
      }
    }
  };

#ifdef USE_OMP
  #pragma omp parallel
  {
    int thread_id = omp_get_thread_num();
    efficient_twoe(thread_id);
  }
#else
  efficient_twoe(0);
#endif
  for(int i = 0; i < this->controls_->nthreads; i++) AX += G[i];
//RealMatrix Tmp = 0.5*(AX + AX.transpose());
//AX = 0.25*Tmp; // Can't consolidate where this comes from?
//if(doRHFFock) AX = 0.25*AX; // Can't consolidate where this comes from?
  AX = AX*0.5; // Gaussian nonsense
  AX = AX*0.5; // werid factor that comes from A + AT
  if(doRHFFock) AX = AX*0.5; // E ~ 0.5*G
  finish = std::chrono::high_resolution_clock::now();
  if(doRHFFock) this->PTD = finish - start;
   
}
template<>
void AOIntegrals::twoEContractN4(bool doRHFFock, const RealMatrix &X, RealMatrix &AX) {
  this->fileio_->out << "Contracting with in-core two-electron integrals" << endl;
  if(!this->haveAOTwoE) this->computeAOTwoE();
  RealTensor2d XTensor(X.rows(),X.cols());
  RealTensor2d AXTensor(AX.rows(),AX.cols());
  for(auto i = 0; i < X.size(); i++) XTensor.storage()[i] = X.data()[i];
  AXTensor.fill(0.0);

  double fact = -1.0;
  if(doRHFFock) fact = -0.5;

  enum{i,j,k,l}; 
  contract(1.0,*this->aoERI_,{i,j,k,l},XTensor,{l,k},0.0,AXTensor,{i,j});
  contract(fact,*this->aoERI_,{i,l,k,j},XTensor,{l,k},1.0,AXTensor,{i,j});

  for(auto i = 0; i < AX.size(); i++) AX.data()[i] = AXTensor.storage()[i];
//AX = AX*0.5; // Gaussian nonsense
//AX = AX*0.5; // werid factor that comes from A + AT
  if(doRHFFock) AX = AX*0.5; // E ~ 0.5*G
}
#endif

}
