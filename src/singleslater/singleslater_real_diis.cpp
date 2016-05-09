/* 
 *  The Chronus Quantum (ChronusQ) software package is high-performace 
 *  computational chemistry software with a strong emphasis on explicitly 
 *  time-dependent and post-SCF quantum mechanical methods.
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
#include <singleslater.h>
using ChronusQ::FileIO;
using ChronusQ::SingleSlater;

namespace ChronusQ{
/*
template<>
void SingleSlater<double>::CDIIS(){
  int N = this->lenCoeff_;
  int NRHS = 1;
  int INFO = -1;
  int NB = this->nBasis_ * this->nTCS_; 
  int NBSq = NB*NB;
  char NORM = 'O';
  double RCOND;

//cout << "BEFORE ALLOC" << endl;
//this->memManager_->printSummary(cout);
  RealMap B(this->memManager_->malloc<double>(N*N),N,N);
  double * coef   = this->memManager_->malloc<double>(N);
  int    * iPiv   = this->memManager_->malloc<int>(N);
  int    * iWORK_ = this->memManager_->malloc<int>(N);
//cout << "After ALLOC" << endl;
//this->memManager_->printSummary(cout);

  for(auto j = 0; j < (N-1); j++)
  for(auto k = 0; k <= j          ; k++){
    RealMap EJA(this->ErrorAlphaMem_ + (j%(N-1))*NBSq,NB,NB);
    RealMap EKA(this->ErrorAlphaMem_ + (k%(N-1))*NBSq,NB,NB);
    B(j,k) = -EJA.frobInner(EKA);
    if(!this->isClosedShell && this->Ref_ != TCS){
      RealMap EJB(this->ErrorBetaMem_ + (j%(N-1))*NBSq,NB,NB);
      RealMap EKB(this->ErrorBetaMem_ + (k%(N-1))*NBSq,NB,NB);
      B(j,k) += -EJB.frobInner(EKB);
    }
    B(k,j) = B(j,k);
  }

  for (auto l=0;l<N-1;l++){
     B(N-1,l)=-1.0;
     B(l,N-1)=-1.0;
  }

  B(N-1,N-1)=0;
  for(auto k = 0; k < N;k++) coef[k] = 0.0; 
  coef[N-1]=-1.0;

  double ANORM = B.lpNorm<1>();
  //std::vector<int> iWORK_(N);
  dgesv_(&N,&NRHS,B.data(),&N,iPiv,coef,&N,&INFO);
  //dgecon_(&NORM,&N,B.data(),&N,&ANORM,&RCOND,this->WORK_,&iWORK_[0],&INFO);
  dgecon_(&NORM,&N,B.data(),&N,&ANORM,&RCOND,this->WORK_,iWORK_,&INFO);


  if(std::abs(RCOND) > std::numeric_limits<double>::epsilon()) {
    this->fockA_->setZero();
    if(!this->isClosedShell && this->nTCS_ == 1) 
      this->fockB_->setZero();

    for(auto j = 0; j < N-1; j++) {
      RealMap FA(this->FADIIS_ + (j%(N-1))*NBSq,NB,NB);
      (*this->fockA_) += coef[j]*FA;
      if(!this->isClosedShell && this->Ref_ != TCS) {
        RealMap FB(this->FBDIIS_ + (j%(N-1))*NBSq,NB,NB);
        (*this->fockB_) += coef[j]*FB;
      }
    }
  }
  this->memManager_->free(B.data(),N*N);
  this->memManager_->free(coef,N);
  this->memManager_->free(iPiv,N);
  this->memManager_->free(iWORK_,N);
//cout << "After FREE" << endl;
//this->memManager_->printSummary(cout);

} // CDIIS
*/

template<>
void SingleSlater<double>::CpyFock(int iter){
  auto NTCSxNBASIS = this->nTCS_ * this->nBasis_;
  auto NSQ = NTCSxNBASIS  * NTCSxNBASIS;
  std::memcpy(this->FADIIS_+(iter % (this->nDIISExtrap_-1)) * NSQ,
    this->fockA_->data(),NSQ * sizeof(double));

  if(!this->isClosedShell && this->Ref_ != TCS)
    std::memcpy(this->FBDIIS_ + (iter % (this->nDIISExtrap_-1)) * NSQ,
                this->fockB_->data(),NSQ * sizeof(double));
} // CpyFock

template<>
void SingleSlater<double>::GenDComm(int iter){
  RealMap ErrA(
    this->ErrorAlphaMem_ + (iter % (this->nDIISExtrap_-1)) * this->lenF_,
    this->nTCS_*this->nBasis_,this->nTCS_*this->nBasis_
  );

  ErrA = (*this->fockA_) * (*this->onePDMA_) * (*this->aointegrals_->overlap_);
  ErrA -= (*this->aointegrals_->overlap_) * (*this->onePDMA_) * (*this->fockA_);
  if(!this->isClosedShell && this->Ref_ != TCS){
    RealMap ErrB(
      this->ErrorBetaMem_ + (iter % (this->nDIISExtrap_-1)) * this->lenF_,
      this->nBasis_,this->nBasis_
    );

    ErrB = (*this->fockB_) * (*this->onePDMB_) * (*this->aointegrals_->overlap_);
    ErrB -= (*this->aointegrals_->overlap_) * (*this->onePDMB_) * (*this->fockB_);
  }
} // GenDComm

template<>
void SingleSlater<double>::genDComm2(int iter) {
  auto NTCSxNBASIS = this->nTCS_ * this->nBasis_;
  auto NSQ = NTCSxNBASIS  * NTCSxNBASIS;
  RealMap ErrA(
    this->ErrorAlphaMem_ + (iter % (this->nDIISExtrap_-1)) * NSQ,
    this->nTCS_*this->nBasis_,this->nTCS_*this->nBasis_
  );

  if(this->nTCS_ == 1 && this->isClosedShell) {
    (*this->NBSqScratch_) = (*this->onePDMA_) * (*this->aointegrals_->overlap_);
    ErrA = (*this->fockA_) * (*this->NBSqScratch_);
    ErrA -= this->NBSqScratch_->adjoint() * (*this->fockA_);
  } else {
    RealMap ErrB(
      this->ErrorBetaMem_ + (iter % (this->nDIISExtrap_-1)) * NSQ,
      this->nTCS_*this->nBasis_,this->nTCS_*this->nBasis_
    );
    (*this->NBSqScratch_) = (*this->onePDMB_) * (*this->aointegrals_->overlap_);
    ErrB = (*this->fockB_) * (*this->NBSqScratch_);
    ErrB -= this->NBSqScratch_->adjoint() * (*this->fockB_);
  };
};

}// Namespace ChronusQ
