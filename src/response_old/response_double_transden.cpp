/* 
 *  This file is part of the Chronus Quantum (ChronusQ) software package
 *  
 *  Copyright (C) 2014-2016 Li Research Group (University of Washington)
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

namespace ChronusQ {

/**
 *  Places the vectorized contents of T into the virtual-occupied block
 *  of TMOA and TMOB
 */
template<>
void Response<double>::placeVirOcc(RealVecMap &T, RealMatrix &TMOA,
  RealMatrix &TMOB) {

  bool doBeta = this->iPart_ != SPIN_ADAPTED && 
                this->Ref_ != SingleSlater<double>::TCS &&
                this->iClass_ != PPPA;

  Eigen::Map<RealMatrix>  TExpandedA(T.data(),this->nVA_,this->nOA_);
  Eigen::Map<RealMatrix>  TExpandedB(T.data(),0,0);
  if(doBeta)
    new (&TExpandedB) Eigen::Map<RealMatrix>(
      T.data()+this->nOAVA_,this->nVB_,this->nOB_
    );
  
  auto iOffA = this->nOA_;
  auto iOffB = this->nOB_;

  TMOA.block(iOffA,0,this->nVA_,this->nOA_) = TExpandedA;
  if(doBeta)
    TMOB.block(iOffB,0,this->nVA_,this->nOA_) = TExpandedB;

}; // placeVirOcc

template<>
void Response<double>::placeOccVir(RealVecMap &T, RealMatrix &TMOA,
  RealMatrix &TMOB) {

  bool doBeta = this->iPart_ != SPIN_ADAPTED && 
                this->Ref_ != SingleSlater<double>::TCS &&
                this->iClass_ != PPPA;

  Eigen::Map<RealMatrix>  TExpandedA(T.data(),this->nVA_,this->nOA_);
  Eigen::Map<RealMatrix>  TExpandedB(T.data(),0,0);
  if(doBeta)
    new (&TExpandedB) Eigen::Map<RealMatrix>(
      T.data()+this->nOAVA_,this->nVB_,this->nOB_
    );
  
  auto iOffA = this->nOA_;
  auto iOffB = this->nOB_;

  TMOA.block(0,iOffA,this->nOA_,this->nVA_) = TExpandedA.adjoint();
  if(doBeta)
    TMOB.block(0,iOffB,this->nOA_,this->nVA_) = TExpandedB.adjoint();

}; // placeVirOcc

template<>
void Response<double>::retrieveVirOcc(RealVecMap &T, RealMatrix &TMOA,
  RealMatrix &TMOB) {

  bool doBeta = this->iPart_ != SPIN_ADAPTED && 
                this->Ref_ != SingleSlater<double>::TCS &&
                this->iClass_ != PPPA;

  auto iOffA = this->nOA_;
  auto iOffB = this->nOB_;

  RealMatrix TExpandedA, TExpandedB;
  TExpandedA = TMOA.block(iOffA,0,this->nVA_,this->nOA_);
  if(doBeta) TExpandedB = TMOB.block(iOffB,0,this->nVB_,this->nOB_);

  RealVecMap TLinA(TExpandedA.data(),0);
  RealVecMap TLinB(TExpandedA.data(),0);

  new (&TLinA) RealVecMap(TExpandedA.data(),this->nOAVA_);
  if(doBeta) new (&TLinB) RealVecMap(TExpandedB.data(),this->nOBVB_);

  T.block(0,0,this->nOAVA_,1) = TLinA;
  if(doBeta) T.block(this->nOAVA_,0,this->nOBVB_,1) = TLinB;

} // retrieveVirOcc

template<>
void Response<double>::retrieveOccVir(RealVecMap &T, RealMatrix &TMOA,
  RealMatrix &TMOB) {

  bool doBeta = this->iPart_ != SPIN_ADAPTED && 
                this->Ref_ != SingleSlater<double>::TCS &&
                this->iClass_ != PPPA;

  auto iOffA = this->nOA_;
  auto iOffB = this->nOB_;

  RealMatrix TExpandedA, TExpandedB;
  TExpandedA = TMOA.block(0,iOffA,this->nOA_,this->nVA_);
  if(doBeta) TExpandedB = TMOB.block(0,iOffB,this->nOB_,this->nVB_);

  TExpandedA.transposeInPlace();
  if(doBeta) TExpandedB.transposeInPlace();

  RealVecMap TLinA(TExpandedA.data(),0);
  RealVecMap TLinB(TExpandedA.data(),0);

  new (&TLinA) RealVecMap(TExpandedA.data(),this->nOAVA_);
  if(doBeta) new (&TLinB) RealVecMap(TExpandedB.data(),this->nOBVB_);

  // Wierd sign that I can't consolidate
  T.block(0,0,this->nOAVA_,1) = -TLinA;
  if(doBeta) T.block(this->nOAVA_,0,this->nOBVB_,1) = -TLinB;

} // retrieveOccVir

template<>
void Response<double>::formAOTransDenFOPPA(RealVecMap &T, RealMatrix &TAOA,
  RealMatrix &TAOB) {
  RealMatrix TMOA,TMOB;
  bool doBeta = this->iPart_ != SPIN_ADAPTED && 
                this->Ref_ != SingleSlater<double>::TCS &&
                this->iClass_ != PPPA;

  TMOA = RealMatrix(this->nBasis_,this->nBasis_);
  if(doBeta)
    TMOB = RealMatrix(this->nBasis_,this->nBasis_);

  this->placeVirOcc(T,TMOA,TMOB);
  if(!this->doTDA_){
    RealVecMap Y(
      T.data()+this->nSingleDim_/2,this->nSingleDim_/2
    );

    this->placeOccVir(Y,TMOA,TMOB);
  }

  TAOA = (*this->singleSlater_->moA()) * TMOA * 
         this->singleSlater_->moA()->adjoint();
  if(doBeta){
    if(this->singleSlater_->isClosedShell)
      TAOB = (*this->singleSlater_->moA()) * TMOB * 
             this->singleSlater_->moA()->adjoint();
    else
      TAOB = (*this->singleSlater_->moB()) * TMOB * 
             this->singleSlater_->moB()->adjoint();
  }
}; //formAOTDenFOPPA

template<>
void Response<double>::formMOTransDenFOPPA(RealVecMap &T, RealMatrix &TAOA,
  RealMatrix &TAOB) {
  RealMatrix TMOA,TMOB;
  bool doBeta = this->iPart_ != SPIN_ADAPTED && 
                this->Ref_ != SingleSlater<double>::TCS &&
                this->iClass_ != PPPA;

  TMOA = RealMatrix(this->nBasis_,this->nBasis_);
  if(doBeta)
    TMOB = RealMatrix(this->nBasis_,this->nBasis_);

  TMOA = this->singleSlater_->moA()->adjoint() * TAOA * 
         (*this->singleSlater_->moA());
  if(doBeta){
    if(this->singleSlater_->isClosedShell)
      TMOB = this->singleSlater_->moA()->adjoint() * TAOB * 
             (*this->singleSlater_->moA());
    else
      TMOB = this->singleSlater_->moB()->adjoint() * TAOB * 
             (*this->singleSlater_->moB());
  }

  this->retrieveVirOcc(T,TMOA,TMOB);
  if(!this->doTDA_){
    RealVecMap Y(
      T.data()+this->nSingleDim_/2,this->nSingleDim_/2
    );

    this->retrieveOccVir(Y,TMOA,TMOB);
  }

}; //formMOTDenFOPPA

template<>
void Response<double>::formAOTransDenPPRPA(RealVecMap &T, RealMatrix &TAOA,
  RealMatrix &TAOB) {

  for(auto mu = 0; mu < this->nTCS_*this->nBasis_; mu++)
  for(auto nu = 0; nu < this->nTCS_*this->nBasis_; nu++){
    if(this->currentMat_ == AA_PPRPA || this->currentMat_ == AAA_PPTDA ||
       this->currentMat_ == PPRPA_TRIPLETS) {
      for(auto a = 0, ab = 0; a < this->nVA_; a++      )
      for(auto b = 0        ; b < a        ;  b++, ab++){
        TAOA(mu,nu) += (*this->singleSlater_->moA())(mu,a) *
                       (*this->singleSlater_->moA())(nu,b) *
                       T(ab);
      }
    } 
    // FIXME: Need to generalize for UHF
    else if(this->currentMat_ == AB_PPRPA || this->currentMat_ == AAB_PPTDA){
      for(auto a = 0, ab = 0; a < this->nVA_; a++      )
      for(auto b = 0        ; b < this->nVB_; b++, ab++){
        TAOA(mu,nu) += (*this->singleSlater_->moA())(mu,a) *
                       (*this->singleSlater_->moA())(nu,b) *
                       T(ab);
      }
    }
  }
}; //formAOTransDenPPRPA

template<>
void Response<double>::formMOTransDenPPRPA(RealVecMap &T, RealMatrix &TMOA,
  RealMatrix &TMOB) {

}; //formMOTransDenPPRPA

}; // namespace ChronusQ
