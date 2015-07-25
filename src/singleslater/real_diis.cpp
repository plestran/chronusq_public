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
#include <singleslater.h>
using ChronusQ::Controls;
using ChronusQ::FileIO;
using ChronusQ::SingleSlater;

namespace ChronusQ{
template<>
void SingleSlater<double>::CDIIS(int N, double *EA, double *FADIIS, double *EB, double *FBDIIS){
  RealMatrix B(N,N);
  double *coef = new double[N];
  int    *iPiv = new int[N];
  int    NRHS = 1, INFO = -1;
  int NBSq = this->nBasis_*this->nBasis_;
  for(auto j = 0; j < (N-1); j++)
  for(auto k = 0; k <= j          ; k++){
    RealMap EJA(EA + (j%(N-1))*NBSq,this->nBasis_,this->nBasis_);
    if((k==0) && this->controls_->printLevel > 4) 
      prettyPrint(this->fileio_->out,EJA,"Error "+std::to_string(j));
    RealMap EKA(EA + (k%(N-1))*NBSq,this->nBasis_,this->nBasis_);
    B(j,k) = -EJA.frobInner(EKA);
    if(this->Ref_ != RHF){
      RealMap EJB(EB + (j%(N-1))*NBSq,this->nBasis_,this->nBasis_);
      RealMap EKB(EB + (k%(N-1))*NBSq,this->nBasis_,this->nBasis_);
      B(j,k) = -EJB.frobInner(EKB);
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
  dgesv_(&N,&NRHS,B.data(),&N,iPiv,coef,&N,&INFO);
  this->fockA_->setZero();
  if(this->Ref_ != RHF) this->fockB_->setZero();
  for(auto j = 0; j < N-1; j++) {
    RealMap FA(FADIIS + (j%(N-1))*NBSq,this->nBasis_,this->nBasis_);
    *this->fockA_ += coef[j]*FA;
    if(this->Ref_ != RHF) {
      RealMap FB(FBDIIS + (j%(N-1))*NBSq,this->nBasis_,this->nBasis_);
      *this->fockB_ += coef[j]*FB;
    }
  }
  if(this->controls_->printLevel > 4) prettyPrint(this->fileio_->out,*this->fockA_,"Total Fock");
  delete [] coef;
  delete [] iPiv;

}
}// Namespace ChronusQ
