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
void AOIntegrals::twoEContractDirect(bool RHF, bool doFock, const RealMatrix &XAlpha, RealMatrix &AXAlpha,
                                     const RealMatrix &XBeta, RealMatrix &AXBeta) {
  this->fileio_->out << "Contracting Directly with two-electron integrals" << endl;

  if(!this->haveSchwartz) this->computeSchwartz();
  if(!this->basisSet_->haveMap) this->basisSet_->makeMap(this->molecule_); 
  AXAlpha.setZero();
  AXBeta.setZero();
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
                    if(RHF && doFock) {
                      G[0][thread_id](bf1,bf2) += XAlpha(bf3,bf4)*v;
                      G[0][thread_id](bf3,bf4) += XAlpha(bf1,bf2)*v;
                      G[0][thread_id](bf2,bf1) += XAlpha(bf4,bf3)*v;
                      G[0][thread_id](bf4,bf3) += XAlpha(bf2,bf1)*v;
                    } else if(doFock) {
                      G[0][thread_id](bf1,bf2) += (XAlpha(bf3,bf4)+XBeta(bf3,bf4))*v;
                      G[0][thread_id](bf3,bf4) += (XAlpha(bf1,bf2)+XBeta(bf1,bf2))*v;
                      G[0][thread_id](bf2,bf1) += (XAlpha(bf4,bf3)+XBeta(bf4,bf3))*v;
                      G[0][thread_id](bf4,bf3) += (XAlpha(bf2,bf1)+XBeta(bf2,bf1))*v;
                      G[1][thread_id](bf1,bf2) += (XAlpha(bf3,bf4)+XBeta(bf3,bf4))*v;
                      G[1][thread_id](bf3,bf4) += (XAlpha(bf1,bf2)+XBeta(bf1,bf2))*v;
                      G[1][thread_id](bf2,bf1) += (XAlpha(bf4,bf3)+XBeta(bf4,bf3))*v;
                      G[1][thread_id](bf4,bf3) += (XAlpha(bf2,bf1)+XBeta(bf2,bf1))*v;
                    } else {
                      G[0][thread_id](bf1,bf2) += (XAlpha(bf3,bf4)+XBeta(bf3,bf4))*v;
                      G[0][thread_id](bf3,bf4) += (XAlpha(bf1,bf2)+XBeta(bf1,bf2))*v;
                      G[0][thread_id](bf2,bf1) += (XAlpha(bf4,bf3)+XBeta(bf4,bf3))*v;
                      G[0][thread_id](bf4,bf3) += (XAlpha(bf2,bf1)+XBeta(bf2,bf1))*v;
                      G[1][thread_id](bf1,bf2) += (XAlpha(bf3,bf4)+XBeta(bf3,bf4))*v;
                      G[1][thread_id](bf3,bf4) += (XAlpha(bf1,bf2)+XBeta(bf1,bf2))*v;
                      G[1][thread_id](bf2,bf1) += (XAlpha(bf4,bf3)+XBeta(bf4,bf3))*v;
                      G[1][thread_id](bf4,bf3) += (XAlpha(bf2,bf1)+XBeta(bf2,bf1))*v;

                      G[0][thread_id](bf1,bf2) += (XAlpha(bf4,bf3)+XBeta(bf4,bf3))*v;
                      G[0][thread_id](bf3,bf4) += (XAlpha(bf2,bf1)+XBeta(bf2,bf1))*v;
                      G[0][thread_id](bf2,bf1) += (XAlpha(bf3,bf4)+XBeta(bf3,bf4))*v;
                      G[0][thread_id](bf4,bf3) += (XAlpha(bf1,bf2)+XBeta(bf1,bf2))*v;
                      G[1][thread_id](bf1,bf2) += (XAlpha(bf4,bf3)+XBeta(bf4,bf3))*v;
                      G[1][thread_id](bf3,bf4) += (XAlpha(bf2,bf1)+XBeta(bf2,bf1))*v;
                      G[1][thread_id](bf2,bf1) += (XAlpha(bf3,bf4)+XBeta(bf3,bf4))*v;
                      G[1][thread_id](bf4,bf3) += (XAlpha(bf1,bf2)+XBeta(bf1,bf2))*v;
                    }
 
                    // Exchange
                    if(RHF && doFock) {
                      G[0][thread_id](bf1,bf3) -= 0.25*XAlpha(bf2,bf4)*v;
                      G[0][thread_id](bf4,bf2) -= 0.25*XAlpha(bf1,bf3)*v;
                      G[0][thread_id](bf1,bf4) -= 0.25*XAlpha(bf2,bf3)*v;
                      G[0][thread_id](bf3,bf2) -= 0.25*XAlpha(bf1,bf4)*v;

                      G[0][thread_id](bf3,bf1) -= 0.25*XAlpha(bf4,bf2)*v;
                      G[0][thread_id](bf2,bf4) -= 0.25*XAlpha(bf3,bf1)*v;
                      G[0][thread_id](bf4,bf1) -= 0.25*XAlpha(bf3,bf2)*v;
                      G[0][thread_id](bf2,bf3) -= 0.25*XAlpha(bf4,bf1)*v;
                    } else if(doFock) {
                      G[0][thread_id](bf1,bf3) -= 0.5*XAlpha(bf2,bf4)*v;
                      G[0][thread_id](bf4,bf2) -= 0.5*XAlpha(bf1,bf3)*v;
                      G[0][thread_id](bf1,bf4) -= 0.5*XAlpha(bf2,bf3)*v;
                      G[0][thread_id](bf3,bf2) -= 0.5*XAlpha(bf1,bf4)*v;

                      G[0][thread_id](bf3,bf1) -= 0.5*XAlpha(bf4,bf2)*v;
                      G[0][thread_id](bf2,bf4) -= 0.5*XAlpha(bf3,bf1)*v;
                      G[0][thread_id](bf4,bf1) -= 0.5*XAlpha(bf3,bf2)*v;
                      G[0][thread_id](bf2,bf3) -= 0.5*XAlpha(bf4,bf1)*v;

                      G[1][thread_id](bf1,bf3) -= 0.5*XBeta(bf2,bf4)*v;
                      G[1][thread_id](bf4,bf2) -= 0.5*XBeta(bf1,bf3)*v;
                      G[1][thread_id](bf1,bf4) -= 0.5*XBeta(bf2,bf3)*v;
                      G[1][thread_id](bf3,bf2) -= 0.5*XBeta(bf1,bf4)*v;

                      G[1][thread_id](bf3,bf1) -= 0.5*XBeta(bf4,bf2)*v;
                      G[1][thread_id](bf2,bf4) -= 0.5*XBeta(bf3,bf1)*v;
                      G[1][thread_id](bf4,bf1) -= 0.5*XBeta(bf3,bf2)*v;
                      G[1][thread_id](bf2,bf3) -= 0.5*XBeta(bf4,bf1)*v;
                    } else {
                      G[0][thread_id](bf1,bf3) -= XAlpha(bf2,bf4)*v;
                      G[0][thread_id](bf4,bf2) -= XAlpha(bf1,bf3)*v;
                      G[0][thread_id](bf1,bf4) -= XAlpha(bf2,bf3)*v;
                      G[0][thread_id](bf3,bf2) -= XAlpha(bf1,bf4)*v;

                      G[0][thread_id](bf3,bf1) -= XAlpha(bf4,bf2)*v;
                      G[0][thread_id](bf2,bf4) -= XAlpha(bf3,bf1)*v;
                      G[0][thread_id](bf4,bf1) -= XAlpha(bf3,bf2)*v;
                      G[0][thread_id](bf2,bf3) -= XAlpha(bf4,bf1)*v;

                      G[1][thread_id](bf1,bf3) -= XBeta(bf2,bf4)*v;
                      G[1][thread_id](bf4,bf2) -= XBeta(bf1,bf3)*v;
                      G[1][thread_id](bf1,bf4) -= XBeta(bf2,bf3)*v;
                      G[1][thread_id](bf3,bf2) -= XBeta(bf1,bf4)*v;

                      G[1][thread_id](bf3,bf1) -= XBeta(bf4,bf2)*v;
                      G[1][thread_id](bf2,bf4) -= XBeta(bf3,bf1)*v;
                      G[1][thread_id](bf4,bf1) -= XBeta(bf3,bf2)*v;
                      G[1][thread_id](bf2,bf3) -= XBeta(bf4,bf1)*v;
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
#endif

}
