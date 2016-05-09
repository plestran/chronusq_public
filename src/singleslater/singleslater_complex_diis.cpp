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
void SingleSlater<dcomplex>::CDIIS(){
  int N = this->lenCoeff_;
  ComplexMatrix B(N,N);
  dcomplex *coef = new dcomplex[N];
  int    *iPiv = new int[N];
  int    NRHS = 1, INFO = -1;
  int NBSq = this->nBasis_*this->nBasis_*this->nTCS_*this->nTCS_;
  for(auto j = 0; j < (N-1); j++)
  for(auto k = 0; k <= j          ; k++){
    ComplexMap EJA(this->ErrorAlphaMem_ + (j%(N-1))*NBSq,this->nTCS_*this->nBasis_,this->nTCS_*this->nBasis_);
    ComplexMap EKA(this->ErrorAlphaMem_ + (k%(N-1))*NBSq,this->nTCS_*this->nBasis_,this->nTCS_*this->nBasis_);
    B(j,k) = -EJA.frobInner(EKA);
    if(!this->isClosedShell && this->Ref_ != TCS){
      ComplexMap EJB(this->ErrorBetaMem_ + (j%(N-1))*NBSq,this->nBasis_,this->nBasis_);
      ComplexMap EKB(this->ErrorBetaMem_ + (k%(N-1))*NBSq,this->nBasis_,this->nBasis_);
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

  ComplexVecMap COEFF(coef,N);
  VectorXcd     RHS(COEFF);
  COEFF = B.fullPivLu().solve(RHS);

  this->fockA_->setZero();
  if(!this->isClosedShell && this->Ref_ != TCS) this->fockB_->setZero();
  for(auto j = 0; j < N-1; j++) {
    ComplexMap FA(this->FADIIS_ + (j%(N-1))*NBSq,this->nTCS_*this->nBasis_,this->nTCS_*this->nBasis_);
    *this->fockA_ += coef[j]*FA;
    if(!this->isClosedShell && this->Ref_ != TCS) {
      ComplexMap FB(this->FBDIIS_ + (j%(N-1))*NBSq,this->nBasis_,this->nBasis_);
      *this->fockB_ += coef[j]*FB;
    }
  }
  delete [] coef;
  delete [] iPiv;

} // CDIIS
*/

template<>
void SingleSlater<dcomplex>::CpyFock(int iter){
  auto NTCSxNBASIS = this->nTCS_ * this->nBasis_;
  auto NSQ = NTCSxNBASIS  * NTCSxNBASIS;
  std::memcpy(this->FADIIS_+(iter % (this->lenCoeff_-1)) * NSQ,this->fockA_->data(),
              NSQ * sizeof(dcomplex));
  if(!this->isClosedShell && this->Ref_ != TCS)
    std::memcpy(this->FBDIIS_ + (iter % (this->lenCoeff_-1)) * NSQ,
                this->fockB_->data(),NSQ * sizeof(dcomplex));
} // CpyFock

template<>
void SingleSlater<dcomplex>::GenDComm(int iter){
  ComplexMap ErrA(this->ErrorAlphaMem_ + (iter % (this->lenCoeff_-1)) * this->lenF_,
               this->nTCS_*this->nBasis_,this->nTCS_*this->nBasis_);
  ErrA = (*this->fockA_) * (*this->onePDMA_) * (*this->aointegrals_->overlap_);
  ErrA -= (*this->aointegrals_->overlap_) * (*this->onePDMA_) * (*this->fockA_);
  if(!this->isClosedShell && this->Ref_ != TCS){
    ComplexMap ErrB(this->ErrorBetaMem_ + (iter % (this->lenCoeff_-1)) * this->lenF_,
                 this->nBasis_,this->nBasis_);

    ErrB = (*this->fockB_) * (*this->onePDMB_) * (*this->aointegrals_->overlap_);
    ErrB -= (*this->aointegrals_->overlap_) * (*this->onePDMB_) * (*this->fockB_);
  }
} // GenDComm

template<>
void SingleSlater<dcomplex>::genDComm2(int iter) {
  auto NTCSxNBASIS = this->nTCS_ * this->nBasis_;
  auto NSQ = NTCSxNBASIS  * NTCSxNBASIS;
  ComplexMap ErrA(
    this->ErrorAlphaMem_ + (iter % (this->nDIISExtrap_-1)) * NSQ,
    this->nTCS_*this->nBasis_,this->nTCS_*this->nBasis_
  );

  if(this->nTCS_ == 1 && this->isClosedShell) {
    this->NBSqScratch_->real() = 
      this->onePDMA_->real() * (*this->aointegrals_->overlap_);
    this->NBSqScratch_->imag() = 
      this->onePDMA_->imag() * (*this->aointegrals_->overlap_);

    ErrA = (*this->fockA_) * (*this->NBSqScratch_);
    ErrA -= this->NBSqScratch_->adjoint() * (*this->fockA_);
  } else {
    ComplexMap ErrB(
      this->ErrorBetaMem_ + (iter % (this->nDIISExtrap_-1)) * NSQ,
      this->nTCS_*this->nBasis_,this->nTCS_*this->nBasis_
    );

    this->NBSqScratch_->real() = 
      this->onePDMB_->real() * (*this->aointegrals_->overlap_);
    this->NBSqScratch_->imag() = 
      this->onePDMB_->imag() * (*this->aointegrals_->overlap_);

    ErrB = (*this->fockB_) * (*this->NBSqScratch_);
    ErrB -= this->NBSqScratch_->adjoint() * (*this->fockB_);
  };
};

}// Namespace ChronusQ
