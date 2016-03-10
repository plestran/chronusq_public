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
using ChronusQ::Controls;
using ChronusQ::FileIO;
using ChronusQ::SingleSlater;

namespace ChronusQ{
template<>
void SingleSlater<double>::CDIIS(){
  int N = this->lenCoeff_;
  RealMatrix B(N,N);
  double *coef = new double[N];
  int    *iPiv = new int[N];
  int    NRHS = 1, INFO = -1;
  int NBSq = this->nBasis_*this->nBasis_*this->nTCS_*this->nTCS_;
  for(auto j = 0; j < (N-1); j++)
  for(auto k = 0; k <= j          ; k++){
    RealMap EJA(this->ErrorAlphaMem_ + (j%(N-1))*NBSq,this->nTCS_*this->nBasis_,this->nTCS_*this->nBasis_);
    RealMap EKA(this->ErrorAlphaMem_ + (k%(N-1))*NBSq,this->nTCS_*this->nBasis_,this->nTCS_*this->nBasis_);
    B(j,k) = -EJA.frobInner(EKA);
    if(!this->isClosedShell && this->Ref_ != TCS){
      RealMap EJB(this->ErrorBetaMem_ + (j%(N-1))*NBSq,this->nBasis_,this->nBasis_);
      RealMap EKB(this->ErrorBetaMem_ + (k%(N-1))*NBSq,this->nBasis_,this->nBasis_);
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
/*
  prettyPrint(this->fileio_->out,B,"B");
  this->fileio_->out << "COEFF" << endl;
  for(auto k = 0; k < N;k++) this->fileio_->out << coef[k] << endl;
  this->fileio_->out << endl; 
*/

//dgesv_(&N,&NRHS,B.data(),&N,iPiv,coef,&N,&INFO);
  char NORM = 'O';
  double ANORM = B.lpNorm<1>();
//RealVecMap COEFF(coef,N);
//VectorXd   RHS(COEFF);
//COEFF = B.fullPivLu().solve(RHS);

//dgetrf_(&N,&N,B.data(),&N,iPiv,&INFO);
//std::vector<int> iWORK_(N);
//dgecon_(&NORM,&N,B.data(),&N,&ANORM,&RCOND,this->WORK_,&iWORK_[0],&INFO);

  char TRANS = 'N';
  dgels_(&TRANS,&N,&N,&NRHS,B.data(),&N,coef,&N,this->WORK_,&this->LWORK_,&INFO);

  /*
  double RCOND = -1.0;
  int Rank;
  double* S    = new double[N];
  dgelss_(&N,&N,&NRHS,B.data(),&N,coef,&N,S,&RCOND,&Rank,this->WORK_,
    &this->LWORK_,&INFO);
  delete [] S;
  */


//if(std::abs(RCOND) > 1e-10) {
    this->fockA_->setZero();
    if(!this->isClosedShell && this->Ref_ != TCS) this->fockB_->setZero();
    for(auto j = 0; j < N-1; j++) {
      RealMap FA(this->FADIIS_ + (j%(N-1))*NBSq,this->nTCS_*this->nBasis_,this->nTCS_*this->nBasis_);
      *this->fockA_ += coef[j]*FA;
      if(!this->isClosedShell && this->Ref_ != TCS) {
        RealMap FB(this->FBDIIS_ + (j%(N-1))*NBSq,this->nBasis_,this->nBasis_);
        *this->fockB_ += coef[j]*FB;
      }
    }
//}
  delete [] coef;
  delete [] iPiv;

} // CDIIS

template<>
void SingleSlater<double>::CpyFock(int iter){
  std::memcpy(this->FADIIS_+(iter % (this->lenCoeff_-1)) * this->lenF_,this->fockA_->data(),
              this->lenF_ * sizeof(double));
  if(!this->isClosedShell && this->Ref_ != TCS)
    std::memcpy(this->FBDIIS_ + (iter % (this->lenCoeff_-1)) * this->lenF_,
                this->fockB_->data(),this->lenF_ * sizeof(double));
} // CpyFock

template<>
void SingleSlater<double>::GenDComm(int iter){
  RealMap ErrA(this->ErrorAlphaMem_ + (iter % (this->lenCoeff_-1)) * this->lenF_,
               this->nTCS_*this->nBasis_,this->nTCS_*this->nBasis_);
/*
  if(this->Ref_ == TCS){
    RealMap GenOverlap(this->SMem_,this->nTCS_*this->nBasis_,this->nTCS_*this->nBasis_);
    ErrA = (*this->fockA_) * (*this->densityA_) * (GenOverlap);
    ErrA -= (GenOverlap) * (*this->densityA_) * (*this->fockA_);
  } else {
    ErrA = (*this->fockA_) * (*this->densityA_) * (*this->aointegrals_->overlap_);
    ErrA -= (*this->aointegrals_->overlap_) * (*this->densityA_) * (*this->fockA_);
  }
*/
  ErrA = (*this->fockA_) * (*this->densityA_) * (*this->aointegrals_->overlap_);
  ErrA -= (*this->aointegrals_->overlap_) * (*this->densityA_) * (*this->fockA_);
  if(!this->isClosedShell && this->Ref_ != TCS){
    RealMap ErrB(this->ErrorBetaMem_ + (iter % (this->lenCoeff_-1)) * this->lenF_,
                 this->nBasis_,this->nBasis_);

    ErrB = (*this->fockB_) * (*this->densityB_) * (*this->aointegrals_->overlap_);
    ErrB -= (*this->aointegrals_->overlap_) * (*this->densityB_) * (*this->fockB_);
  }
} // GenDComm

}// Namespace ChronusQ
