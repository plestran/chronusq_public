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

template<typename T>
void SingleSlater<T>::CDIIS(){
  int N = this->lenCoeff_;
  int NRHS = 1;
  int INFO = -1;
  int NB = this->nBasis_ * this->nTCS_; 
  int NBSq = NB*NB;
  char NORM = 'O';
  double RCOND;

  TMap B(this->memManager_->template malloc<T>(N*N),N,N);
  T   * coef   = this->memManager_->template malloc<T>(N);
  int * iPiv   = this->memManager_->template malloc<int>(N);
  int * iWORK_ = this->memManager_->template malloc<int>(N);

  for(auto j = 0; j < (N-1); j++)
  for(auto k = 0; k <= j          ; k++){
    TMap EJA(this->ErrorAlphaMem_ + (j%(N-1))*NBSq,NB,NB);
    TMap EKA(this->ErrorAlphaMem_ + (k%(N-1))*NBSq,NB,NB);
    B(j,k) = -EJA.frobInner(EKA);
    if(!this->isClosedShell && this->Ref_ != TCS){
      TMap EJB(this->ErrorBetaMem_ + (j%(N-1))*NBSq,NB,NB);
      TMap EKB(this->ErrorBetaMem_ + (k%(N-1))*NBSq,NB,NB);
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

  double ANORM = B.template lpNorm<1>();

  if(typeid(T).hash_code() == typeid(dcomplex).hash_code()){
    zgesv_(&N,&NRHS,reinterpret_cast<dcomplex*>(B.data()),&N,iPiv,
        reinterpret_cast<dcomplex*>(coef),&N,&INFO);
    zgecon_(&NORM,&N,reinterpret_cast<dcomplex*>(B.data()),&N,&ANORM,&RCOND,
        reinterpret_cast<dcomplex*>(this->WORK_),this->RWORK_,&INFO);
  } else {
    dgesv_(&N,&NRHS,reinterpret_cast<double*>(B.data()),&N,iPiv,
        reinterpret_cast<double*>(coef),&N,&INFO);
    dgecon_(&NORM,&N,reinterpret_cast<double*>(B.data()),&N,&ANORM,&RCOND,
        reinterpret_cast<double*>(this->WORK_),iWORK_,&INFO);
  }


  if(std::abs(RCOND) > std::numeric_limits<double>::epsilon()) {
    this->fockA_->setZero();
    if(!this->isClosedShell && this->nTCS_ == 1) 
      this->fockB_->setZero();

    for(auto j = 0; j < N-1; j++) {
      TMap FA(this->FADIIS_ + (j%(N-1))*NBSq,NB,NB);
      (*this->fockA_) += coef[j]*FA;
      if(!this->isClosedShell && this->Ref_ != TCS) {
        TMap FB(this->FBDIIS_ + (j%(N-1))*NBSq,NB,NB);
        (*this->fockB_) += coef[j]*FB;
      }
    }
  }
  this->memManager_->free(B.data(),N*N);
  this->memManager_->free(coef,N);
  this->memManager_->free(iPiv,N);
  this->memManager_->free(iWORK_,N);

} // CDIIS
