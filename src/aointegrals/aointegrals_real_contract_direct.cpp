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
  template<>
  void AOIntegrals::twoEContractDirect(bool RHF, bool doFock, bool do24, bool doTCS, 
    const RealMatrix &XAlpha, RealMatrix &AXAlpha, const RealMatrix &XBeta, 
    RealMatrix &AXBeta) {

    int nTCS = 1;
    if(doTCS) nTCS = 2;
//  this->fileio_->out << "Contracting Directly with two-electron integrals" << endl;
    if(!this->haveSchwartz) this->computeSchwartz();
    if(!this->basisSet_->haveMapSh2Bf) this->basisSet_->makeMapSh2Bf(nTCS); 
    AXAlpha.setZero();
    if(!RHF && !doTCS) AXBeta.setZero();
    int nRHF;
    if(RHF || doTCS) nRHF = 1;
    else    nRHF = 2;
    std::vector<std::vector<RealMatrix>> G(nRHF,std::vector<RealMatrix>
      (omp_get_max_threads(),RealMatrix::Zero(nTCS*this->nBasis_,nTCS*this->nBasis_)));
  
    RealMatrix XTotal;
    if(!RHF && !doTCS) {
      XTotal = XAlpha + XBeta;
      if(!doFock) XTotal = 0.5*XTotal + 0.5*(XAlpha.adjoint() + XBeta.adjoint());
    }
  
    std::vector<coulombEngine> engines(omp_get_max_threads());
    engines[0] = coulombEngine(this->basisSet_->maxPrim(),this->basisSet_->maxL(),0);
    engines[0].set_precision(std::numeric_limits<double>::epsilon());
  
    for(int i=1; i<omp_get_max_threads(); i++) engines[i] = engines[0];
  
    auto start = std::chrono::high_resolution_clock::now();
    this->basisSet_->computeShBlkNorm(!RHF && !doTCS,nTCS,&XAlpha,&XBeta);
    auto finish = std::chrono::high_resolution_clock::now();
    if(doFock) this->DenShBlkD = finish - start;
    int ijkl = 0;
    start = std::chrono::high_resolution_clock::now();
  
    auto efficient_twoe = [&] (int thread_id) {
      coulombEngine &engine = engines[thread_id];
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
              if(s1234 % omp_get_max_threads() != thread_id) continue;
              int bf4_s = this->basisSet_->mapSh2Bf(s4);
              int n4    = this->basisSet_->shells(s4).size();
        
              // Schwartz and Density screening
              double shMax;
              if(RHF || doTCS){
                shMax = std::max((*this->basisSet_->shBlkNormAlpha)(s1,s4),
                        std::max((*this->basisSet_->shBlkNormAlpha)(s2,s4),
                        std::max((*this->basisSet_->shBlkNormAlpha)(s3,s4),
                        std::max((*this->basisSet_->shBlkNormAlpha)(s1,s3),
                        std::max((*this->basisSet_->shBlkNormAlpha)(s2,s3),
                                 (*this->basisSet_->shBlkNormAlpha)(s1,s2)))))) * 
                        (*this->schwartz_)(s1,s2) * (*this->schwartz_)(s3,s4);
              } else {
                shMax = std::max((*this->basisSet_->shBlkNormAlpha)(s1,s4),
                        std::max((*this->basisSet_->shBlkNormAlpha)(s2,s4),
                        std::max((*this->basisSet_->shBlkNormAlpha)(s3,s4),
                        std::max((*this->basisSet_->shBlkNormAlpha)(s1,s3),
                        std::max((*this->basisSet_->shBlkNormAlpha)(s2,s3),
                        std::max((*this->basisSet_->shBlkNormAlpha)(s1,s2), 
                        std::max((*this->basisSet_->shBlkNormBeta)(s1,s4),
                        std::max((*this->basisSet_->shBlkNormBeta)(s2,s4),
                        std::max((*this->basisSet_->shBlkNormBeta)(s3,s4),
                        std::max((*this->basisSet_->shBlkNormBeta)(s1,s3),
                        std::max((*this->basisSet_->shBlkNormBeta)(s2,s3),
                                 (*this->basisSet_->shBlkNormBeta)(s1,s2)))))))))))) * 
                        (*this->schwartz_)(s1,s2) * (*this->schwartz_)(s3,s4);
              }
  
              if(shMax < this->controls_->thresholdSchawrtz ) continue;
   
              const double* buff = engine.compute(
                this->basisSet_->shells(s1),
                this->basisSet_->shells(s2),
                this->basisSet_->shells(s3),
                this->basisSet_->shells(s4));
        
              double s12_deg = (s1 == s2) ? 1.0 : 2.0;
              double s34_deg = (s3 == s4) ? 1.0 : 2.0;
              double s12_34_deg = (s1 == s3) ? (s2 == s4 ? 1.0 : 2.0) : 2.0;
              double s1234_deg = s12_deg * s34_deg * s12_34_deg;
  
              if(RHF && doFock) 
                this->Restricted34Contract(G[0][thread_id],XAlpha,n1,n2,n3,n4,
                  bf1_s,bf2_s,bf3_s,bf4_s,buff,s1234_deg);
              else if(!do24 && !doTCS)
                this->UnRestricted34Contract(G[0][thread_id],XAlpha,G[1][thread_id],
                  XBeta,XTotal,n1,n2,n3,n4,bf1_s,bf2_s,bf3_s,bf4_s,buff,s1234_deg);
              else if(doTCS && !do24)
                this->Spinor34Contract(G[0][thread_id],XAlpha,n1,n2,n3,n4,
                  bf1_s,bf2_s,bf3_s,bf4_s,buff,s1234_deg);
              else if(!doTCS && do24)
                this->General24CouContract(G[0][thread_id],XAlpha,n1,n2,n3,n4,
                  bf1_s,bf2_s,bf3_s,bf4_s,buff,s1234_deg);
              else if(doTCS && do24)
                this->Spinor24CouContract(G[0][thread_id],XAlpha,n1,n2,n3,n4,
                  bf1_s,bf2_s,bf3_s,bf4_s,buff,s1234_deg);
            }
          }
        }
      }
    };
  
  
  #ifdef _OPENMP
    #pragma omp parallel
    {
      int thread_id = omp_get_thread_num();
      efficient_twoe(thread_id);
    }
  #else
    efficient_twoe(0);
  #endif
    for(int i = 0; i < omp_get_max_threads(); i++) AXAlpha += G[0][i];
    if(!RHF && !doTCS) for(int i = 0; i < omp_get_max_threads(); i++) AXBeta += G[1][i];
    AXAlpha = AXAlpha*0.5; // werid factor that comes from A + AT
    AXAlpha = AXAlpha*0.5;
    if(do24) AXAlpha *= 0.5;
    if(!RHF && !doTCS){
      AXBeta = AXBeta*0.5; // werid factor that comes from A + AT
      AXBeta = AXBeta*0.5;
      if(do24) AXBeta *= 0.5;
    }
    finish = std::chrono::high_resolution_clock::now();
    if(doFock) this->PTD = finish - start;
     
  }
  template<>
  void AOIntegrals::multTwoEContractDirect(int nVec, bool RHF, bool doFock, bool do24, 
    bool doTCS, const std::vector<RealMatrix> &XAlpha, std::vector<RealMatrix> &AXAlpha, 
    const std::vector<RealMatrix> &XBeta, std::vector<RealMatrix> &AXBeta) {

//  this->fileio_->out << "Contracting Directly with two-electron integrals" << endl;
  
    int nTCS = 1;
    if(doTCS) nTCS = 2;
    if(!this->haveSchwartz) this->computeSchwartz();
    if(!this->basisSet_->haveMapSh2Bf) this->basisSet_->makeMapSh2Bf(nTCS); 
    for(auto i = 0; i < nVec; i++) AXAlpha[i].setZero();
    if(!RHF && !doTCS)
      for(auto i = 0; i < nVec; i++) AXBeta[i].setZero();
    int nRHF;
    if(RHF || doTCS) nRHF = 1;
    else    nRHF = 2;
  std::vector<std::vector<std::vector<RealMatrix>>> G(nRHF,
    std::vector<std::vector<RealMatrix>>(nVec,
      std::vector<RealMatrix>(omp_get_max_threads(),
        RealMatrix::Zero(nTCS*this->nBasis_,nTCS*this->nBasis_))));
  
    std::vector<RealMatrix> XTotal;
    if(!RHF && !doTCS){
      for(auto iX = 0; iX < nVec; iX++){
        XTotal.push_back(XAlpha[iX] + XBeta[iX]);
        if(!doFock)
          XTotal[iX] = 0.5*XTotal[iX] + 0.5*(XAlpha[iX].adjoint() + XBeta[iX].adjoint());
      }
    }
  
    std::vector<coulombEngine> engines(omp_get_max_threads());
    engines[0] = coulombEngine(this->basisSet_->maxPrim(),this->basisSet_->maxL(),0);
    engines[0].set_precision(std::numeric_limits<double>::epsilon());
  
    for(int i=1; i<omp_get_max_threads(); i++) engines[i] = engines[0];
  
    std::vector<RealMatrix> alphaShBlk;
    std::vector<RealMatrix> betaShBlk;
    auto start = std::chrono::high_resolution_clock::now();
    for(auto i = 0; i < nVec; i++){
      this->basisSet_->computeShBlkNorm(!RHF && !doTCS,nTCS,&XAlpha[i],&XBeta[i]);
      alphaShBlk.push_back(*this->basisSet_->shBlkNormAlpha);
      if(!RHF && !doTCS)
        betaShBlk.push_back(*this->basisSet_->shBlkNormBeta);
    }
    auto finish = std::chrono::high_resolution_clock::now();
    if(doFock) this->DenShBlkD = finish - start;
    int ijkl = 0;
    start = std::chrono::high_resolution_clock::now();
  
    auto efficient_twoe = [&] (int thread_id) {
      coulombEngine &engine = engines[thread_id];
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
              if(s1234 % omp_get_max_threads() != thread_id) continue;
              int bf4_s = this->basisSet_->mapSh2Bf(s4);
              int n4    = this->basisSet_->shells(s4).size();
        
              // Schwartz and Density screening
              double shMax;
              for(auto k = 0; k < nVec; k++){
                double tmpMax = 0;
                if(RHF || doTCS){
                  tmpMax = std::max(alphaShBlk[k](s1,s4), std::max(alphaShBlk[k](s2,s4),
                           std::max(alphaShBlk[k](s3,s4), std::max(alphaShBlk[k](s1,s3),
                           std::max(alphaShBlk[k](s2,s3), alphaShBlk[k](s1,s2)))))) * 
                           (*this->schwartz_)(s1,s2) * (*this->schwartz_)(s3,s4);
                } else {
                  tmpMax = std::max(alphaShBlk[k](s1,s4), std::max(alphaShBlk[k](s2,s4),
                           std::max(alphaShBlk[k](s3,s4), std::max(alphaShBlk[k](s1,s3),
                           std::max(alphaShBlk[k](s2,s3), std::max(alphaShBlk[k](s1,s2), 
                           std::max(betaShBlk[k](s1,s4), std::max(betaShBlk[k](s2,s4),
                           std::max(betaShBlk[k](s3,s4), std::max(betaShBlk[k](s1,s3),
                           std::max(betaShBlk[k](s2,s3), betaShBlk[k](s1,s2)))))))))))) * 
                           (*this->schwartz_)(s1,s2) * (*this->schwartz_)(s3,s4);
                }
                if(k == 0 || tmpMax > shMax) shMax = tmpMax;
              }
  
              if(shMax < this->controls_->thresholdSchawrtz ) continue;
   
              const double* buff = engine.compute(
                this->basisSet_->shells(s1),
                this->basisSet_->shells(s2),
                this->basisSet_->shells(s3),
                this->basisSet_->shells(s4));
        
              double s12_deg = (s1 == s2) ? 1.0 : 2.0;
              double s34_deg = (s3 == s4) ? 1.0 : 2.0;
              double s12_34_deg = (s1 == s3) ? (s2 == s4 ? 1.0 : 2.0) : 2.0;
              double s1234_deg = s12_deg * s34_deg * s12_34_deg;
  
              for(auto iX = 0; iX < nVec; iX++){
                if(RHF && doFock) 
                  this->Restricted34Contract(G[0][iX][thread_id],XAlpha[iX],n1,n2,n3,n4,
                    bf1_s,bf2_s,bf3_s,bf4_s,buff,s1234_deg);
                else if(!do24 && !doTCS)
                  this->UnRestricted34Contract(G[0][iX][thread_id],XAlpha[iX],
                    G[1][iX][thread_id],XBeta[iX],XTotal[iX],n1,n2,n3,n4,bf1_s,bf2_s,bf3_s,
                    bf4_s,buff,s1234_deg);
                else if(doTCS && !do24)
                  this->Spinor34Contract(G[0][iX][thread_id],XAlpha[iX],n1,n2,n3,n4,
                    bf1_s,bf2_s,bf3_s,bf4_s,buff,s1234_deg);
                else if(!doTCS && do24)
                  this->General24CouContract(G[0][iX][thread_id],XAlpha[iX],n1,n2,n3,n4,
                    bf1_s,bf2_s,bf3_s,bf4_s,buff,s1234_deg);
                else if(doTCS && do24)
                  this->Spinor24CouContract(G[0][iX][thread_id],XAlpha[iX],n1,n2,n3,n4,
                    bf1_s,bf2_s,bf3_s,bf4_s,buff,s1234_deg);
              } // Loop iX
            } // Loop s4
          } // Loop s3
        } // Loop s2
      } // Loop s1
    };// efficient_twoe 
  
  
  #ifdef _OPENMP
    #pragma omp parallel
    {
      int thread_id = omp_get_thread_num();
      efficient_twoe(thread_id);
    }
  #else
    efficient_twoe(0);
  #endif
    for(int k = 0; k < nVec; k++)
    for(int i = 0; i < omp_get_max_threads(); i++) 
      AXAlpha[k] += G[0][k][i];
    if(!RHF && !doTCS)
      for(int k = 0; k < nVec; k++)
      for(int i = 0; i < omp_get_max_threads(); i++) 
        AXBeta[k] += G[1][k][i];
  
    double fact = 0.25;
    if(do24) fact *= 0.5;
    for(auto k = 0; k < nVec; k++) AXAlpha[k] *= fact;
    if(!RHF && !doTCS) for(auto k = 0; k < nVec; k++) AXBeta[k] *= fact;
    finish = std::chrono::high_resolution_clock::now();
    if(doFock) this->PTD = finish - start;
     
  }
}; // namespace ChronusQ
