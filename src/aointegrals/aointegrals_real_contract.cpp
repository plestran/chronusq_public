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
/**
 *  Given a shell quartet block of ERIs (row-major), the following contraction is performed
 *
 *  G(μ,v) = (μ v | λ σ) X(σ,λ) - 0.5 (μ σ | λ v) X(σ,λ)
 *
 *  This assumes that G is actually an arbitrary spin block of G and that X is actually
 *  X-Total, i.e. X = X-Alpha + X-Beta (hence the 0.5 on the exchange part)
 */
template<>
void AOIntegrals::Restricted34HerContract(RealMatrix &G, const RealMatrix &X, int n1, int n2, int n3, int n4, 
                     int bf1_s, int bf2_s, int bf3_s, int bf4_s, const double* buff, double deg){
  for(int i = 0, ijkl = 0 ; i < n1; ++i) {
    int bf1 = bf1_s + i;
    for(int j = 0; j < n2; ++j) {
      int bf2 = bf2_s + j;
      for(int k = 0; k < n3; ++k) {
        int bf3 = bf3_s + k;
        for(int l = 0; l < n4; ++l, ++ijkl) {
          int bf4 = bf4_s + l;
          double v = buff[ijkl]*deg;

          // Coulomb
          G(bf1,bf2) += X(bf4,bf3)*v;
          G(bf3,bf4) += X(bf2,bf1)*v;
          G(bf2,bf1) += X(bf3,bf4)*v;
          G(bf4,bf3) += X(bf1,bf2)*v;

          // Exchange
          G(bf1,bf3) -= 0.25*X(bf2,bf4)*v;
          G(bf2,bf4) -= 0.25*X(bf1,bf3)*v;
          G(bf1,bf4) -= 0.25*X(bf2,bf3)*v;
          G(bf2,bf3) -= 0.25*X(bf1,bf4)*v;

          G(bf3,bf1) -= 0.25*X(bf4,bf2)*v;
          G(bf4,bf2) -= 0.25*X(bf3,bf1)*v;
          G(bf4,bf1) -= 0.25*X(bf3,bf2)*v;
          G(bf3,bf2) -= 0.25*X(bf4,bf1)*v;
        }
      }
    }
  }
} // Restricted34HerContract 
template<>
void AOIntegrals::UnRestricted34HerContract(RealMatrix &GAlpha, const RealMatrix &XAlpha, RealMatrix &GBeta, const RealMatrix &XBeta, int n1, int n2, int n3, int n4, 
                     int bf1_s, int bf2_s, int bf3_s, int bf4_s, const double* buff, double deg){
  for(int i = 0, ijkl = 0 ; i < n1; ++i) {
    int bf1 = bf1_s + i;
    for(int j = 0; j < n2; ++j) {
      int bf2 = bf2_s + j;
      for(int k = 0; k < n3; ++k) {
        int bf3 = bf3_s + k;
        for(int l = 0; l < n4; ++l, ++ijkl) {
          int bf4 = bf4_s + l;
          double v = buff[ijkl]*deg;

          // Coulomb
          GAlpha(bf1,bf2) += (XAlpha(bf4,bf3)+XBeta(bf4,bf3))*v;
          GAlpha(bf3,bf4) += (XAlpha(bf2,bf1)+XBeta(bf2,bf1))*v;
          GAlpha(bf2,bf1) += (XAlpha(bf3,bf4)+XBeta(bf3,bf4))*v;
          GAlpha(bf4,bf3) += (XAlpha(bf1,bf2)+XBeta(bf1,bf2))*v;
          GBeta(bf1,bf2)  += (XAlpha(bf4,bf3)+XBeta(bf4,bf3))*v;
          GBeta(bf3,bf4)  += (XAlpha(bf2,bf1)+XBeta(bf2,bf1))*v;
          GBeta(bf2,bf1)  += (XAlpha(bf3,bf4)+XBeta(bf3,bf4))*v;
          GBeta(bf4,bf3)  += (XAlpha(bf1,bf2)+XBeta(bf1,bf2))*v;

          // Exchange
          GAlpha(bf1,bf3) -= 0.5*XAlpha(bf2,bf4)*v;
          GAlpha(bf2,bf4) -= 0.5*XAlpha(bf1,bf3)*v;
          GAlpha(bf1,bf4) -= 0.5*XAlpha(bf2,bf3)*v;
          GAlpha(bf2,bf3) -= 0.5*XAlpha(bf1,bf4)*v;

          GAlpha(bf3,bf1) -= 0.5*XAlpha(bf4,bf2)*v;
          GAlpha(bf4,bf2) -= 0.5*XAlpha(bf3,bf1)*v;
          GAlpha(bf4,bf1) -= 0.5*XAlpha(bf3,bf2)*v;
          GAlpha(bf3,bf2) -= 0.5*XAlpha(bf4,bf1)*v;

          GBeta(bf1,bf3)  -= 0.5*XBeta(bf2,bf4)*v;
          GBeta(bf2,bf4)  -= 0.5*XBeta(bf1,bf3)*v;
          GBeta(bf1,bf4)  -= 0.5*XBeta(bf2,bf3)*v;
          GBeta(bf2,bf3)  -= 0.5*XBeta(bf1,bf4)*v;

          GBeta(bf3,bf1)  -= 0.5*XBeta(bf4,bf2)*v;
          GBeta(bf4,bf2)  -= 0.5*XBeta(bf3,bf1)*v;
          GBeta(bf4,bf1)  -= 0.5*XBeta(bf3,bf2)*v;
          GBeta(bf3,bf2)  -= 0.5*XBeta(bf4,bf1)*v;
        }
      }
    }
  }
} // UnRestricted34HerContract
template<>
void AOIntegrals::General34NonHerContract(RealMatrix &GAlpha, const RealMatrix &XAlpha, RealMatrix &GBeta, const RealMatrix &XBeta, int n1, int n2, int n3, int n4, 
                     int bf1_s, int bf2_s, int bf3_s, int bf4_s, const double* buff, double deg){
  for(int i = 0, ijkl = 0 ; i < n1; ++i) {
    int bf1 = bf1_s + i;
    for(int j = 0; j < n2; ++j) {
      int bf2 = bf2_s + j;
      for(int k = 0; k < n3; ++k) {
        int bf3 = bf3_s + k;
        for(int l = 0; l < n4; ++l, ++ijkl) {
          int bf4 = bf4_s + l;
          double v = buff[ijkl]*deg;

          // Coulomb
          GAlpha(bf1,bf2) += 0.5*(XAlpha(bf4,bf3)+XBeta(bf4,bf3))*v;
          GAlpha(bf3,bf4) += 0.5*(XAlpha(bf2,bf1)+XBeta(bf2,bf1))*v;
          GAlpha(bf2,bf1) += 0.5*(XAlpha(bf3,bf4)+XBeta(bf3,bf4))*v;
          GAlpha(bf4,bf3) += 0.5*(XAlpha(bf1,bf2)+XBeta(bf1,bf2))*v;
          GBeta(bf1,bf2)  += 0.5*(XAlpha(bf4,bf3)+XBeta(bf4,bf3))*v;
          GBeta(bf3,bf4)  += 0.5*(XAlpha(bf2,bf1)+XBeta(bf2,bf1))*v;
          GBeta(bf2,bf1)  += 0.5*(XAlpha(bf3,bf4)+XBeta(bf3,bf4))*v;
          GBeta(bf4,bf3)  += 0.5*(XAlpha(bf1,bf2)+XBeta(bf1,bf2))*v;

          GAlpha(bf1,bf2) += 0.5*(XAlpha(bf3,bf4)+XBeta(bf3,bf4))*v;
          GAlpha(bf3,bf4) += 0.5*(XAlpha(bf1,bf2)+XBeta(bf1,bf2))*v;
          GAlpha(bf2,bf1) += 0.5*(XAlpha(bf4,bf3)+XBeta(bf4,bf3))*v;
          GAlpha(bf4,bf3) += 0.5*(XAlpha(bf2,bf1)+XBeta(bf2,bf1))*v;
          GBeta(bf1,bf2)  += 0.5*(XAlpha(bf3,bf4)+XBeta(bf3,bf4))*v;
          GBeta(bf3,bf4)  += 0.5*(XAlpha(bf1,bf2)+XBeta(bf1,bf2))*v;
          GBeta(bf2,bf1)  += 0.5*(XAlpha(bf4,bf3)+XBeta(bf4,bf3))*v;
          GBeta(bf4,bf3)  += 0.5*(XAlpha(bf2,bf1)+XBeta(bf2,bf1))*v;

          // Exchange
          GAlpha(bf1,bf3) -= 0.5*XAlpha(bf2,bf4)*v;
          GAlpha(bf2,bf4) -= 0.5*XAlpha(bf1,bf3)*v;
          GAlpha(bf1,bf4) -= 0.5*XAlpha(bf2,bf3)*v;
          GAlpha(bf2,bf3) -= 0.5*XAlpha(bf1,bf4)*v;

          GAlpha(bf3,bf1) -= 0.5*XAlpha(bf4,bf2)*v;
          GAlpha(bf4,bf2) -= 0.5*XAlpha(bf3,bf1)*v;
          GAlpha(bf4,bf1) -= 0.5*XAlpha(bf3,bf2)*v;
          GAlpha(bf3,bf2) -= 0.5*XAlpha(bf4,bf1)*v;

          GBeta(bf1,bf3)  -= 0.5*XBeta(bf2,bf4)*v;
          GBeta(bf2,bf4)  -= 0.5*XBeta(bf1,bf3)*v;
          GBeta(bf1,bf4)  -= 0.5*XBeta(bf2,bf3)*v;
          GBeta(bf2,bf3)  -= 0.5*XBeta(bf1,bf4)*v;

          GBeta(bf3,bf1)  -= 0.5*XBeta(bf4,bf2)*v;
          GBeta(bf4,bf2)  -= 0.5*XBeta(bf3,bf1)*v;
          GBeta(bf4,bf1)  -= 0.5*XBeta(bf3,bf2)*v;
          GBeta(bf3,bf2)  -= 0.5*XBeta(bf4,bf1)*v;
        }
      }
    }
  }
}
template<>
void AOIntegrals::twoEContractDirect(bool RHF, bool doFock, const RealMatrix &XAlpha, RealMatrix &AXAlpha,
                                     const RealMatrix &XBeta, RealMatrix &AXBeta) {
  this->fileio_->out << "Contracting Directly with two-electron integrals" << endl;

  if(!this->haveSchwartz) this->computeSchwartz();
  if(!this->basisSet_->haveMap) this->basisSet_->makeMap(this->molecule_); 
  AXAlpha.setZero();
  if(!RHF) AXBeta.setZero();
  int nRHF;
  if(RHF) nRHF = 1;
  else    nRHF = 2;
  std::vector<std::vector<RealMatrix>> G(nRHF,std::vector<RealMatrix>
    (this->controls_->nthreads,RealMatrix::Zero(this->nBasis_,this->nBasis_)));

  std::vector<coulombEngine> engines(this->controls_->nthreads);
  engines[0] = coulombEngine(this->basisSet_->maxPrim,this->basisSet_->maxL,0);
  engines[0].set_precision(std::numeric_limits<double>::epsilon());

  for(int i=1; i<this->controls_->nthreads; i++) engines[i] = engines[0];

  auto start = std::chrono::high_resolution_clock::now();
  this->basisSet_->computeShBlkNorm(RHF,this->molecule_,&XAlpha,&XBeta);
  auto finish = std::chrono::high_resolution_clock::now();
  if(doFock) this->DenShBlkD = finish - start;
  int ijkl = 0;
  start = std::chrono::high_resolution_clock::now();

  auto efficient_twoe = [&] (int thread_id) {
    coulombEngine &engine = engines[thread_id];
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
            double shMax;
            if(RHF){
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
              this->basisSet_->shells_libint[s1],
              this->basisSet_->shells_libint[s2],
              this->basisSet_->shells_libint[s3],
              this->basisSet_->shells_libint[s4]);
      
            double s12_deg = (s1 == s2) ? 1.0 : 2.0;
            double s34_deg = (s3 == s4) ? 1.0 : 2.0;
            double s12_34_deg = (s1 == s3) ? (s2 == s4 ? 1.0 : 2.0) : 2.0;
            double s1234_deg = s12_deg * s34_deg * s12_34_deg;

            if(RHF && doFock) 
              this->Restricted34HerContract(G[0][thread_id],XAlpha,n1,n2,n3,n4,
                bf1_s,bf2_s,bf3_s,bf4_s,buff,s1234_deg);
            else if(doFock)
              this->UnRestricted34HerContract(G[0][thread_id],XAlpha,G[1][thread_id],
                XBeta,n1,n2,n3,n4,bf1_s,bf2_s,bf3_s,bf4_s,buff,s1234_deg);
            else
              this->General34NonHerContract(G[0][thread_id],XAlpha,G[1][thread_id],
                XBeta,n1,n2,n3,n4,bf1_s,bf2_s,bf3_s,bf4_s,buff,s1234_deg);
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
  for(int i = 0; i < this->controls_->nthreads; i++) AXAlpha += G[0][i];
  if(!RHF) for(int i = 0; i < this->controls_->nthreads; i++) AXBeta += G[1][i];
  AXAlpha = AXAlpha*0.5; // werid factor that comes from A + AT
  AXAlpha = AXAlpha*0.5;
  if(!RHF){
    AXBeta = AXBeta*0.5; // werid factor that comes from A + AT
    AXBeta = AXBeta*0.5;
  }
  finish = std::chrono::high_resolution_clock::now();
  if(doFock) this->PTD = finish - start;
   
}
template<>
void AOIntegrals::twoEContractN4(bool RHF, bool doFock, const RealMatrix &XAlpha, RealMatrix &AXAlpha,
                                 const RealMatrix &XBeta, RealMatrix &AXBeta) {
  this->fileio_->out << "Contracting with in-core two-electron integrals" << endl;
  if(!this->haveAOTwoE) this->computeAOTwoE();

  RealTensor2d XAlphaTensor,XBetaTensor;
  RealTensor2d AXAlphaTensor,AXBetaTensor;
  enum{i,j,k,l}; 

  if(doFock)  {
    XAlphaTensor  = RealTensor2d(XAlpha.rows(),XAlpha.cols());
    AXAlphaTensor = RealTensor2d(AXAlpha.rows(),AXAlpha.cols());
    for(auto i = 0; i < XAlpha.size(); i++) XAlphaTensor.storage()[i] = XAlpha.data()[i];
    AXAlphaTensor.fill(0.0);
    if(RHF){
      contract(1.0,*this->aoERI_,{i,j,k,l},XAlphaTensor,{l,k},0.0,AXAlphaTensor,{i,j});
      contract(-0.5,*this->aoERI_,{i,l,k,j},XAlphaTensor,{l,k},1.0,AXAlphaTensor,{i,j});
     
    } else {
      XBetaTensor  = RealTensor2d(XBeta.rows(),XBeta.cols());
      AXBetaTensor = RealTensor2d(AXBeta.rows(),AXBeta.cols());
      for(auto i = 0; i < XBeta.size(); i++) XBetaTensor.storage()[i] = XBeta.data()[i];
      AXBetaTensor.fill(0.0);
      RealTensor2d XTotalTensor(XAlpha.rows(),XBeta.cols());

      XTotalTensor = XAlphaTensor + XBetaTensor;

      contract(1.0,*this->aoERI_,{i,j,k,l},XTotalTensor,{l,k},0.0,AXAlphaTensor,{i,j});
      AXBetaTensor = AXAlphaTensor;
      contract(-1.0,*this->aoERI_,{i,l,k,j},XAlphaTensor,{l,k},1.0,AXAlphaTensor,{i,j});
      contract(-1.0,*this->aoERI_,{i,l,k,j},XBetaTensor,{l,k},1.0,AXBetaTensor,{i,j});
    } 
  } else CErr("General Contraction NYI for in-core integrals");
   for(auto i = 0; i < AXAlpha.size(); i++) AXAlpha.data()[i] = AXAlphaTensor.storage()[i];
   if(!RHF)
     for(auto i = 0; i < AXBeta.size(); i++) AXBeta.data()[i] = AXBetaTensor.storage()[i];
}

template<>
void AOIntegrals::twoEContractDF(bool RHF, bool doFock, const RealMatrix &XAlpha, RealMatrix &AXAlpha,
                                 const RealMatrix &XBeta, RealMatrix &AXBeta) {
  this->fileio_->out << "Contracting with in-core density fitting integrals" << endl;
  if(!this->haveRIS)  this->computeAORIS();
  if(!this->haveRII)  this->computeAORII();
  if(!this->haveTRII) this->transformAORII();

  RealTensor2d XAlphaTensor(XAlpha.rows(),XAlpha.cols());
  RealTensor2d AXAlphaTensor(AXAlpha.rows(),AXAlpha.cols());

  for(auto i = 0; i < XAlpha.size(); i++) XAlphaTensor.storage()[i] = XAlpha.data()[i];
  AXAlphaTensor.fill(0.0);

  RealTensor1d A(this->DFbasisSet_->nBasis());
  RealTensor3d B(this->nBasis_,this->nBasis_,this->DFbasisSet_->nBasis());

  double fact = -1.0;
  if(RHF && doFock) fact = -0.5;

  enum{i,j,k,l,XX,YY}; 
  contract(1.0,*this->aoRII_,{i,j,XX},XAlphaTensor,{i,j},0.0,A,{XX});
  contract(1.0,*this->aoRII_,{k,i,XX},XAlphaTensor,{k,j},0.0,B,{i,j,XX});
  contract(1.0,*this->aoRII_,{i,j,XX},A,{XX},0.0,AXAlphaTensor,{i,j});
  contract(-0.5,*this->aoRII_,{i,k,XX},B,{j,k,XX},1.0,AXAlphaTensor,{i,j});

  for(auto i = 0; i < AXAlpha.size(); i++) AXAlpha.data()[i] = AXAlphaTensor.storage()[i];
}
template<>
void AOIntegrals::multTwoEContractDirect(int nVec, bool RHF, bool doFock, const std::vector<RealMatrix> &XAlpha, std::vector<RealMatrix> &AXAlpha,
                                     const std::vector<RealMatrix> &XBeta, std::vector<RealMatrix> &AXBeta) {
  this->fileio_->out << "Contracting Directly with two-electron integrals" << endl;

  if(!this->haveSchwartz) this->computeSchwartz();
  if(!this->basisSet_->haveMap) this->basisSet_->makeMap(this->molecule_); 
  for(auto i = 0; i < nVec; i++) AXAlpha[i].setZero();
  if(!RHF)
    for(auto i = 0; i < nVec; i++) AXBeta[i].setZero();
  int nRHF;
  if(RHF) nRHF = 1;
  else    nRHF = 2;
std::vector<std::vector<std::vector<RealMatrix>>> G(nRHF,std::vector<std::vector<RealMatrix>>(nVec,
  std::vector<RealMatrix>(this->controls_->nthreads,
  RealMatrix::Zero(this->nBasis_,this->nBasis_))));

  std::vector<coulombEngine> engines(this->controls_->nthreads);
  engines[0] = coulombEngine(this->basisSet_->maxPrim,this->basisSet_->maxL,0);
  engines[0].set_precision(std::numeric_limits<double>::epsilon());

  for(int i=1; i<this->controls_->nthreads; i++) engines[i] = engines[0];

  std::vector<RealMatrix> alphaShBlk;
  std::vector<RealMatrix> betaShBlk;
  auto start = std::chrono::high_resolution_clock::now();
  for(auto i = 0; i < nVec; i++){
    this->basisSet_->computeShBlkNorm(RHF,this->molecule_,&XAlpha[i],&XBeta[i]);
    alphaShBlk.push_back(*this->basisSet_->shBlkNormAlpha);
    if(!RHF)
      betaShBlk.push_back(*this->basisSet_->shBlkNormBeta);
  }
  auto finish = std::chrono::high_resolution_clock::now();
  if(doFock) this->DenShBlkD = finish - start;
  int ijkl = 0;
  start = std::chrono::high_resolution_clock::now();

  auto efficient_twoe = [&] (int thread_id) {
    coulombEngine &engine = engines[thread_id];
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
            double shMax;
            for(auto k = 0; k < nVec; k++){
              double tmpMax = 0;
              if(RHF){
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
              this->basisSet_->shells_libint[s1],
              this->basisSet_->shells_libint[s2],
              this->basisSet_->shells_libint[s3],
              this->basisSet_->shells_libint[s4]);
      
            double s12_deg = (s1 == s2) ? 1.0 : 2.0;
            double s34_deg = (s3 == s4) ? 1.0 : 2.0;
            double s12_34_deg = (s1 == s3) ? (s2 == s4 ? 1.0 : 2.0) : 2.0;
            double s1234_deg = s12_deg * s34_deg * s12_34_deg;

            for(auto iX = 0; iX < nVec; iX++){
              if(RHF && doFock) 
                this->Restricted34HerContract(G[0][iX][thread_id],XAlpha[iX],n1,n2,n3,n4,
                  bf1_s,bf2_s,bf3_s,bf4_s,buff,s1234_deg);
              else if(doFock)
                this->UnRestricted34HerContract(G[0][iX][thread_id],XAlpha[iX],
                  G[1][iX][thread_id],XBeta[iX],n1,n2,n3,n4,bf1_s,bf2_s,bf3_s,bf4_s,buff,
                  s1234_deg);
              else
                this->General34NonHerContract(G[0][iX][thread_id],XAlpha[iX],
                  G[1][iX][thread_id],XBeta[iX],n1,n2,n3,n4,bf1_s,bf2_s,bf3_s,bf4_s,buff,
                  s1234_deg);
            } // Loop iX
          } // Loop s4
        } // Loop s3
      } // Loop s2
    } // Loop s1
  };// efficient_twoe 


#ifdef USE_OMP
  #pragma omp parallel
  {
    int thread_id = omp_get_thread_num();
    efficient_twoe(thread_id);
  }
#else
  efficient_twoe(0);
#endif
  for(int k = 0; k < nVec; k++)
  for(int i = 0; i < this->controls_->nthreads; i++) 
    AXAlpha[k] += G[0][k][i];
  if(!RHF)
    for(int k = 0; k < nVec; k++)
    for(int i = 0; i < this->controls_->nthreads; i++) 
      AXBeta[k] += G[1][k][i];

  for(auto k = 0; k < nVec; k++) AXAlpha[k] *= 0.25;
  if(!RHF) for(auto k = 0; k < nVec; k++) AXBeta[k] *= 0.25;
  finish = std::chrono::high_resolution_clock::now();
  if(doFock) this->PTD = finish - start;
   
}
#endif

}
