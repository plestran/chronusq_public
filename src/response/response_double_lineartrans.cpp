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
#include <response.h>

namespace ChronusQ{
/*
 *  Forms sigma (and possibly rho) for the linear transfomation of E^(2)
 *  (and possibly S^(2) ) onto trial vectors (or any vector in general)
 *
 *  Adapted from Helgaker, et al. JCP 113, 8908 (2000)
 *
 *  DBWY (2015)
 */
template<>
void Response<double>::linearTransFOPPA(RealCMMap &VR, RealCMMap &VL,
  RealCMMap &SR, RealCMMap &SL, RealCMMap &RR, RealCMMap &RL){

  auto NTCSxNBASIS = this->nTCS_ * this->nBasis_;
  auto nVec = VR.cols();
  if(!this->doTDA_) nVec *= 2;

  std::vector<RealCMMatrix> CommA, CommB, GCommA, GCommB;

  // Storage of commutator
  CommA =
    std::vector<RealCMMatrix>(nVec,RealMatrix::Zero(NTCSxNBASIS,NTCSxNBASIS));

  // Storage of G[commutators]
  GCommA =
    std::vector<RealCMMatrix>(nVec,RealMatrix::Zero(NTCSxNBASIS,NTCSxNBASIS));

  if(this->Ref_ != SingleSlater<double>::TCS && this->iPart_ != SPIN_ADAPTED){
    // Storage of commutator
    CommB =
      std::vector<RealCMMatrix>(nVec,RealMatrix::Zero(NTCSxNBASIS,NTCSxNBASIS));
 
    // Storage of G[commutators]
    GCommB =
      std::vector<RealCMMatrix>(nVec,RealMatrix::Zero(NTCSxNBASIS,NTCSxNBASIS));

  }

  // Storage for S*P
  RealMatrix SDA,SDB;

  // Scaling factor for RHF
  double fact = 1.0;
  if(this->Ref_ == SingleSlater<double>::RHF) fact = 0.5;

  // Offset X -> Y
  int iOff = this->nSingleDim_/2; 


  // Compute S*P FIXME: This can be computed once and stored on disk
  SDA = (*this->singleSlater_->aointegrals()->overlap_) * 
        (*this->singleSlater_->densityA());

  if(!this->singleSlater_->isClosedShell && 
     this->Ref_ != SingleSlater<double>::TCS)
    SDB = (*this->singleSlater_->aointegrals()->overlap_) * 
          (*this->singleSlater_->densityB());

  RealCMMatrix XAAO(NTCSxNBASIS,NTCSxNBASIS);
  RealCMMatrix XBAO(NTCSxNBASIS,NTCSxNBASIS);
  RealCMMatrix SigAOA(NTCSxNBASIS,NTCSxNBASIS);
  RealCMMatrix SigAOB(NTCSxNBASIS,NTCSxNBASIS);
  for (auto idx = 0; idx < nVec; idx++){
    /*
     *  XAO(s) = C(s) * XMO(s) * C(s)**H
     *
     *  XAO(s) - s-spin block of X in the AO basis
     *  XMO(s) - s-spin block of X in the MO basis
     *  C(s)   - s-spin block of the MO coefficients
     *  H      - Adjoint
     */ 

    cout << "HERE " << idx << endl;
    RealVecMap X(VR.data(),0);
    if(idx < VR.cols())
     new (&X) RealVecMap(VR.data()+idx*this->nSingleDim_,this->nSingleDim_);
    else
     new (&X) RealVecMap(
       VL.data()+(idx-VR.cols())*this->nSingleDim_,this->nSingleDim_
     );


    this->formAOTransDen(X,XAAO,XBAO);
    CommA[idx] =  fact * XAAO * SDA; 
    CommA[idx] += fact * SDA.adjoint() * XAAO;
    if(this->iPart_ != SPIN_ADAPTED) {
      if(this->Ref_ == SingleSlater<double>::RHF){ 
        CommB[idx] =  fact * XBAO * SDA;
        CommB[idx] += fact * SDA.adjoint() * XBAO;
      } else if(this->Ref_ != SingleSlater<double>::TCS) {
        CommB[idx] =  fact * XBAO * SDB;
        CommB[idx] += fact * SDB.adjoint() * XBAO;
      }
    }

  }; // idx loop

  if( this->singleSlater_->aointegrals()->integralAlgorithm == 
      AOIntegrals::DIRECT 
      && this->nTCS_ != 2)
    this->singleSlater_->aointegrals()->multTwoEContractDirect(nVec,
      false,false,false,false,(this->nTCS_==2),CommA,GCommA,CommB,GCommB);
  else if(this->singleSlater_->aointegrals()->integralAlgorithm == 
          AOIntegrals::INCORE)
    for(auto idx = 0; idx < nVec; idx++)
      this->singleSlater_->aointegrals()->twoEContractN4(false,false,true,
        false,(this->nTCS_==2),CommA[idx],GCommA[idx],CommB[idx],GCommB[idx]);
  else
    CErr("Integral Contraction logic for SDR is not defined",
      this->fileio_->out);

  for(auto idx = 0; idx < nVec; idx++){
    SigAOA =  
      (*this->singleSlater_->fockA()) * CommA[idx] * 
      (*this->singleSlater_->aointegrals()->overlap_);
    SigAOA -= 
      (*this->singleSlater_->aointegrals()->overlap_) * CommA[idx] * 
      (*this->singleSlater_->fockA());
    SigAOA += fact * GCommA[idx] * SDA.adjoint();
    SigAOA -= fact * SDA * GCommA[idx];

    if(this->iPart_ != SPIN_ADAPTED) {
      if(this->Ref_ == SingleSlater<double>::RHF) {
        SigAOB =  
          (*this->singleSlater_->fockA()) * CommB[idx] * 
          (*this->singleSlater_->aointegrals()->overlap_);
        SigAOB -= 
          (*this->singleSlater_->aointegrals()->overlap_) * CommB[idx] * 
          (*this->singleSlater_->fockA());
        SigAOB += fact * GCommB[idx] * SDA.adjoint();
        SigAOB -= fact * SDA * GCommB[idx];
      } else if(this->Ref_ != SingleSlater<double>::TCS) {
        SigAOB =  
          (*this->singleSlater_->fockB()) * CommB[idx] * 
          (*this->singleSlater_->aointegrals()->overlap_);
        SigAOB -= 
          (*this->singleSlater_->aointegrals()->overlap_) * CommB[idx] * 
          (*this->singleSlater_->fockB());
        SigAOB += fact * GCommB[idx] * SDB.adjoint();
        SigAOB -= fact * SDB * GCommB[idx];
      }
    }

    RealVecMap SVec(SR.data()+idx*this->nSingleDim_,this->nSingleDim_);
    this->formMOTransDen(SVec,SigAOA,SigAOB);
    


  }

}; // linearTransFOPPA

}; // namespace ChronusQ
