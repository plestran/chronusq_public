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
/*
  template<>
  void AOIntegrals::twoEContractDirect(bool RHF, bool doFock, bool do24, bool doTCS, const ComplexMatrix &XAlpha, ComplexMatrix &AXAlpha,
                                   const ComplexMatrix &XBeta, ComplexMatrix &AXBeta) {
  //CErr("No Direct Contraction for Complex Matricies Implemented");
    if(doTCS) CErr("TCS Complex contraction NYI");
    this->fileio_->out << "Contracting Directly with two-electron integrals" << endl;
  
    if(!this->haveSchwartz) this->computeSchwartz();
    if(!this->basisSet_->haveMapSh2Bf) this->basisSet_->makeMapSh2Bf(1); 
    AXAlpha.setZero();
    if(!RHF) AXBeta.setZero();
    int nRHF;
    if(RHF) nRHF = 1;
    else    nRHF = 2;
    std::vector<std::vector<ComplexMatrix>> G(nRHF,std::vector<ComplexMatrix>
      (omp_get_max_threads(),ComplexMatrix::Zero(this->nBasis_,this->nBasis_)));
  
    std::vector<coulombEngine> engines(omp_get_max_threads());
    engines[0] = coulombEngine(this->basisSet_->maxPrim(),this->basisSet_->maxL(),0);
    engines[0].set_precision(std::numeric_limits<double>::epsilon());
  
    for(int i=1; i<omp_get_max_threads(); i++) engines[i] = engines[0];
  
    auto start = std::chrono::high_resolution_clock::now();
    this->basisSet_->computeShBlkNorm(!RHF,1,&XAlpha,&XBeta);
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
                this->basisSet_->shells(s1),
                this->basisSet_->shells(s2),
                this->basisSet_->shells(s3),
                this->basisSet_->shells(s4));
        
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
                        G[0][thread_id](bf1,bf2) += XAlpha(bf4,bf3).real()*v;
                        G[0][thread_id](bf3,bf4) += XAlpha(bf2,bf1).real()*v;
                        G[0][thread_id](bf2,bf1) += XAlpha(bf3,bf4).real()*v;
                        G[0][thread_id](bf4,bf3) += XAlpha(bf1,bf2).real()*v;
                      } else if(doFock) {
                        G[0][thread_id](bf1,bf2) += (XAlpha(bf4,bf3)+XBeta(bf4,bf3)).real()*v;
                        G[0][thread_id](bf3,bf4) += (XAlpha(bf2,bf1)+XBeta(bf2,bf1)).real()*v;
                        G[0][thread_id](bf2,bf1) += (XAlpha(bf3,bf4)+XBeta(bf3,bf4)).real()*v;
                        G[0][thread_id](bf4,bf3) += (XAlpha(bf1,bf2)+XBeta(bf1,bf2)).real()*v;
                        G[1][thread_id](bf1,bf2) += (XAlpha(bf4,bf3)+XBeta(bf4,bf3)).real()*v;
                        G[1][thread_id](bf3,bf4) += (XAlpha(bf2,bf1)+XBeta(bf2,bf1)).real()*v;
                        G[1][thread_id](bf2,bf1) += (XAlpha(bf3,bf4)+XBeta(bf3,bf4)).real()*v;
                        G[1][thread_id](bf4,bf3) += (XAlpha(bf1,bf2)+XBeta(bf1,bf2)).real()*v;
                      } else {
                        G[0][thread_id](bf1,bf2) += 0.5*(XAlpha(bf4,bf3)+XBeta(bf4,bf3))*v;
                        G[0][thread_id](bf3,bf4) += 0.5*(XAlpha(bf2,bf1)+XBeta(bf2,bf1))*v;
                        G[0][thread_id](bf2,bf1) += 0.5*(XAlpha(bf3,bf4)+XBeta(bf3,bf4))*v;
                        G[0][thread_id](bf4,bf3) += 0.5*(XAlpha(bf1,bf2)+XBeta(bf1,bf2))*v;
                        G[1][thread_id](bf1,bf2) += 0.5*(XAlpha(bf4,bf3)+XBeta(bf4,bf3))*v;
                        G[1][thread_id](bf3,bf4) += 0.5*(XAlpha(bf2,bf1)+XBeta(bf2,bf1))*v;
                        G[1][thread_id](bf2,bf1) += 0.5*(XAlpha(bf3,bf4)+XBeta(bf3,bf4))*v;
                        G[1][thread_id](bf4,bf3) += 0.5*(XAlpha(bf1,bf2)+XBeta(bf1,bf2))*v;
  
                        G[0][thread_id](bf1,bf2) += 0.5*(XAlpha(bf3,bf4)+XBeta(bf3,bf4))*v;
                        G[0][thread_id](bf3,bf4) += 0.5*(XAlpha(bf1,bf2)+XBeta(bf1,bf2))*v;
                        G[0][thread_id](bf2,bf1) += 0.5*(XAlpha(bf4,bf3)+XBeta(bf4,bf3))*v;
                        G[0][thread_id](bf4,bf3) += 0.5*(XAlpha(bf2,bf1)+XBeta(bf2,bf1))*v;
                        G[1][thread_id](bf1,bf2) += 0.5*(XAlpha(bf3,bf4)+XBeta(bf3,bf4))*v;
                        G[1][thread_id](bf3,bf4) += 0.5*(XAlpha(bf1,bf2)+XBeta(bf1,bf2))*v;
                        G[1][thread_id](bf2,bf1) += 0.5*(XAlpha(bf4,bf3)+XBeta(bf4,bf3))*v;
                        G[1][thread_id](bf4,bf3) += 0.5*(XAlpha(bf2,bf1)+XBeta(bf2,bf1))*v;
                      }
                      // Exchange
                      if(RHF && doFock) {
                        G[0][thread_id](bf1,bf3) -= 0.25*XAlpha(bf2,bf4)*v;
                        G[0][thread_id](bf2,bf4) -= 0.25*XAlpha(bf1,bf3)*v;
                        G[0][thread_id](bf1,bf4) -= 0.25*XAlpha(bf2,bf3)*v;
                        G[0][thread_id](bf2,bf3) -= 0.25*XAlpha(bf1,bf4)*v;
  
                        G[0][thread_id](bf3,bf1) -= 0.25*XAlpha(bf4,bf2)*v;
                        G[0][thread_id](bf4,bf2) -= 0.25*XAlpha(bf3,bf1)*v;
                        G[0][thread_id](bf4,bf1) -= 0.25*XAlpha(bf3,bf2)*v;
                        G[0][thread_id](bf3,bf2) -= 0.25*XAlpha(bf4,bf1)*v;
                      } else if(doFock) {
                        G[0][thread_id](bf1,bf3) -= 0.5*XAlpha(bf2,bf4)*v;
                        G[0][thread_id](bf2,bf4) -= 0.5*XAlpha(bf1,bf3)*v;
                        G[0][thread_id](bf1,bf4) -= 0.5*XAlpha(bf2,bf3)*v;
                        G[0][thread_id](bf2,bf3) -= 0.5*XAlpha(bf1,bf4)*v;
  
                        G[0][thread_id](bf3,bf1) -= 0.5*XAlpha(bf4,bf2)*v;
                        G[0][thread_id](bf4,bf2) -= 0.5*XAlpha(bf3,bf1)*v;
                        G[0][thread_id](bf4,bf1) -= 0.5*XAlpha(bf3,bf2)*v;
                        G[0][thread_id](bf3,bf2) -= 0.5*XAlpha(bf4,bf1)*v;
  
                        G[1][thread_id](bf1,bf3) -= 0.5*XBeta(bf2,bf4)*v;
                        G[1][thread_id](bf2,bf4) -= 0.5*XBeta(bf1,bf3)*v;
                        G[1][thread_id](bf1,bf4) -= 0.5*XBeta(bf2,bf3)*v;
                        G[1][thread_id](bf2,bf3) -= 0.5*XBeta(bf1,bf4)*v;
  
                        G[1][thread_id](bf3,bf1) -= 0.5*XBeta(bf4,bf2)*v;
                        G[1][thread_id](bf4,bf2) -= 0.5*XBeta(bf3,bf1)*v;
                        G[1][thread_id](bf4,bf1) -= 0.5*XBeta(bf3,bf2)*v;
                        G[1][thread_id](bf3,bf2) -= 0.5*XBeta(bf4,bf1)*v;
                      } else {
                        G[0][thread_id](bf1,bf3) -= 0.5*XAlpha(bf2,bf4)*v;
                        G[0][thread_id](bf2,bf4) -= 0.5*XAlpha(bf1,bf3)*v;
                        G[0][thread_id](bf1,bf4) -= 0.5*XAlpha(bf2,bf3)*v;
                        G[0][thread_id](bf2,bf3) -= 0.5*XAlpha(bf1,bf4)*v;
  
                        G[0][thread_id](bf3,bf1) -= 0.5*XAlpha(bf4,bf2)*v;
                        G[0][thread_id](bf4,bf2) -= 0.5*XAlpha(bf3,bf1)*v;
                        G[0][thread_id](bf4,bf1) -= 0.5*XAlpha(bf3,bf2)*v;
                        G[0][thread_id](bf3,bf2) -= 0.5*XAlpha(bf4,bf1)*v;
  
                        G[1][thread_id](bf1,bf3) -= 0.5*XBeta(bf2,bf4)*v;
                        G[1][thread_id](bf2,bf4) -= 0.5*XBeta(bf1,bf3)*v;
                        G[1][thread_id](bf1,bf4) -= 0.5*XBeta(bf2,bf3)*v;
                        G[1][thread_id](bf2,bf3) -= 0.5*XBeta(bf1,bf4)*v;
  
                        G[1][thread_id](bf3,bf1) -= 0.5*XBeta(bf4,bf2)*v;
                        G[1][thread_id](bf4,bf2) -= 0.5*XBeta(bf3,bf1)*v;
                        G[1][thread_id](bf4,bf1) -= 0.5*XBeta(bf3,bf2)*v;
                        G[1][thread_id](bf3,bf2) -= 0.5*XBeta(bf4,bf1)*v;
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
    if(!RHF) for(int i = 0; i < omp_get_max_threads(); i++) AXBeta += G[1][i];
    AXAlpha = AXAlpha*0.5; // werid factor that comes from A + AT
    AXAlpha = AXAlpha*0.5;
    if(!RHF){
      AXBeta = AXBeta*0.5; // werid factor that comes from A + AT
      AXBeta = AXBeta*0.5;
    }
    finish = std::chrono::high_resolution_clock::now();
    if(doFock) this->PTD = finish - start;
  }
*/
//template<>
//void AOIntegrals::twoEContractN4(bool RHF, bool doFock, bool do24, bool doTCS, const ComplexMatrix &XAlpha, ComplexMatrix &AXAlpha,
//                                 const ComplexMatrix &XBeta, ComplexMatrix &AXBeta) {
//  this->fileio_->out << "Contracting with in-core two-electron integrals" << endl;
//  if(!this->haveAOTwoE) this->computeAOTwoE();
///*
//  RealTensor2d ReXAlphaTensor(XAlpha.rows(),XAlpha.cols());
//  RealTensor2d ImXAlphaTensor(XAlpha.rows(),XAlpha.cols());
//  RealTensor2d ReAXAlphaTensor(AXAlpha.rows(),AXAlpha.cols());
//  RealTensor2d ImAXAlphaTensor(AXAlpha.rows(),AXAlpha.cols());
//  RealTensor2d ReXBetaTensor(XBeta.rows(),XBeta.cols());
//  RealTensor2d ImXBetaTensor(XBeta.rows(),XBeta.cols());
//  RealTensor2d ReAXBetaTensor(AXBeta.rows(),AXBeta.cols());
//  RealTensor2d ImAXBetaTensor(AXBeta.rows(),AXBeta.cols());
//  for(auto i = 0; i < XAlpha.size(); i++) {
//    ReXAlphaTensor.storage()[i] = XAlpha.data()[i].real();
//    ImXAlphaTensor.storage()[i] = XAlpha.data()[i].imag();
//    ReXBetaTensor.storage()[i] = XBeta.data()[i].real();
//    ImXBetaTensor.storage()[i] = XBeta.data()[i].imag();
//  }
//  ReAXAlphaTensor.fill(0.0);
//  ImAXAlphaTensor.fill(0.0);
//  ReAXBetaTensor.fill(0.0);
//  ImAXBetaTensor.fill(0.0);
//
//  double fact = -1.0;
//  if(RHF && doFock) fact = -0.5;
//
//  enum{i,j,k,l}; 
//  contract(1.0,*this->aoERI_,{i,j,k,l},ReXAlphaTensor,{l,k},0.0,ReAXAlphaTensor,{i,j});
//  contract(1.0,*this->aoERI_,{i,j,k,l},ImXAlphaTensor,{l,k},0.0,ImAXAlphaTensor,{i,j});
//  contract(fact,*this->aoERI_,{i,l,k,j},ReXAlphaTensor,{l,k},1.0,ReAXAlphaTensor,{i,j});
//  contract(fact,*this->aoERI_,{i,l,k,j},ImXAlphaTensor,{l,k},1.0,ImAXAlphaTensor,{i,j});
//
//  for(auto i = 0; i < AXAlpha.size(); i++) {
//    AXAlpha.data()[i] = 
//      dcomplex(ReAXAlphaTensor.storage()[i],ImAXAlphaTensor.storage()[i]);
//  }
//*/
//  RealTensor2d ReXAlphaTensor, ReXBetaTensor;
//  RealTensor2d ReAXAlphaTensor,ReAXBetaTensor;
//  RealTensor2d ImXAlphaTensor, ImXBetaTensor;
//  RealTensor2d ImAXAlphaTensor,ImAXBetaTensor;
//  enum{i,j,k,l}; 
//
//  if(doFock)  {
//    ReXAlphaTensor  = RealTensor2d(XAlpha.rows(),XAlpha.cols());
//    ReAXAlphaTensor = RealTensor2d(AXAlpha.rows(),AXAlpha.cols());
//    ImXAlphaTensor  = RealTensor2d(XAlpha.rows(),XAlpha.cols());
//    ImAXAlphaTensor = RealTensor2d(AXAlpha.rows(),AXAlpha.cols());
//    for(auto i = 0; i < XAlpha.size(); i++) {
//      ReXAlphaTensor.storage()[i] = XAlpha.data()[i].real();
//      ImXAlphaTensor.storage()[i] = XAlpha.data()[i].imag();
//    }
//    ReAXAlphaTensor.fill(0.0);
//    ImAXAlphaTensor.fill(0.0);
//    if(RHF){
//      contract(1.0,*this->aoERI_,{i,j,k,l},ReXAlphaTensor,{l,k},0.0,ReAXAlphaTensor,{i,j});
//      contract(1.0,*this->aoERI_,{i,j,k,l},ImXAlphaTensor,{l,k},0.0,ImAXAlphaTensor,{i,j});
//      contract(-0.5,*this->aoERI_,{i,l,k,j},ReXAlphaTensor,{l,k},1.0,ReAXAlphaTensor,{i,j});
//      contract(-0.5,*this->aoERI_,{i,l,k,j},ImXAlphaTensor,{l,k},1.0,ImAXAlphaTensor,{i,j});
//     
//    } else {
//      ReXBetaTensor  = RealTensor2d(XBeta.rows(),XBeta.cols());
//      ReAXBetaTensor = RealTensor2d(AXBeta.rows(),AXBeta.cols());
//      ImXBetaTensor  = RealTensor2d(XBeta.rows(),XBeta.cols());
//      ImAXBetaTensor = RealTensor2d(AXBeta.rows(),AXBeta.cols());
//      for(auto i = 0; i < XBeta.size(); i++) {
//        ReXBetaTensor.storage()[i] = XBeta.data()[i].real();
//        ImXBetaTensor.storage()[i] = XBeta.data()[i].imag();
//      }
//      ReAXBetaTensor.fill(0.0);
//      ImAXBetaTensor.fill(0.0);
//      RealTensor2d ReXTotalTensor(XAlpha.rows(),XBeta.cols());
//      RealTensor2d ImXTotalTensor(XAlpha.rows(),XBeta.cols());
//
//      ReXTotalTensor = ReXAlphaTensor + ReXBetaTensor;
//      ImXTotalTensor = ImXAlphaTensor + ImXBetaTensor;
//
//      contract(1.0,*this->aoERI_,{i,j,k,l},ReXTotalTensor,{l,k},0.0,ReAXAlphaTensor,{i,j});
//      contract(1.0,*this->aoERI_,{i,j,k,l},ImXTotalTensor,{l,k},0.0,ImAXAlphaTensor,{i,j});
//      ReAXBetaTensor = ReAXAlphaTensor;
//      ImAXBetaTensor = ImAXAlphaTensor;
//      contract(-1.0,*this->aoERI_,{i,l,k,j},ReXAlphaTensor,{l,k},1.0,ReAXAlphaTensor,{i,j});
//      contract(-1.0,*this->aoERI_,{i,l,k,j},ImXAlphaTensor,{l,k},1.0,ImAXAlphaTensor,{i,j});
//      contract(-1.0,*this->aoERI_,{i,l,k,j},ReXBetaTensor,{l,k},1.0,ReAXBetaTensor,{i,j});
//      contract(-1.0,*this->aoERI_,{i,l,k,j},ImXBetaTensor,{l,k},1.0,ImAXBetaTensor,{i,j});
//    } 
//  } else CErr("General Contraction NYI for in-core integrals");
//  for(auto i = 0; i < AXAlpha.size(); i++) {
//    AXAlpha.data()[i] = 
//      dcomplex(ReAXAlphaTensor.storage()[i],ImAXAlphaTensor.storage()[i]);
//    if(!RHF) AXBeta.data()[i] = 
//      dcomplex(ReAXBetaTensor.storage()[i],ImAXBetaTensor.storage()[i]);
//  }
//}
//template<>
//void AOIntegrals::twoEContractDF(bool RHF, bool doFock, const ComplexMatrix &XAlpha, ComplexMatrix &AXAlpha,
//                                 const ComplexMatrix &XBeta, ComplexMatrix &AXBeta) {
//  CErr("No Density Fitting Contraction for Complex Matricies Implemented");
//}
  
  }
