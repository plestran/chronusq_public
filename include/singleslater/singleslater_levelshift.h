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


/*
//DBWY and JJRADLER  -- Needs to be rewritten and have conditional for convergence added
//JJRADLER - This function should be called within SCF loop until convergence criterion is met.
//This algorithm follows Mehrotra
template<typename T>
void SingleSlater<T>::levelShift(){
//Construction of Delta matrix
// cout << "In Level Shift" << endl;
  double b = 2.42;	// 2 + .The meaning of the Universe.
  TMatrix deltaA = TMatrix::Zero(this->nBasis_, this->nBasis_);
  for(auto iVir = this->nOccA_; iVir < this->nBasis_; iVir++){
    deltaA(iVir, iVir) = b;
  }
  //Transformation of Delta matrix from MO to AO basis
  TMatrix transC = (*this->moA_).transpose();
  TMatrix dagC = transC.conjugate();
  TMatrix deltaAPrime = dagC * deltaA * (*this->moA_);	//Transformation of delta (MO) to deltaPrime (AO)
  
  //Shift current Fock matrix (F + Delta')
  for(auto iVir = this->nOccA_; iVir < this->nBasis_; iVir++){
    (*this->fockA_)(iVir, iVir) += deltaAPrime(iVir, iVir);	//Sums the diagonal elements of the vir-vir blocks only
  }

  //If CUHF is performed, apply same treatment to Beta MOs
  if(this->Ref_ == CUHF) {
    TMatrix deltaB = TMatrix::Zero(this->nBasis_, this->nBasis_);
    for(auto iVir = this->nOccB_; iVir < this->nBasis_; iVir++){
      deltaB(iVir, iVir) = b;
    }
    //Transformation of Delta matrix from MO to AO basis
    TMatrix transC = (*this->moB_).transpose();
    TMatrix dagC = transC.conjugate();
    TMatrix deltaBPrime = dagC * deltaB * (*this->moB_);	//Transformation of delta (MO) to deltaPrime (AO)
  
    //Shift current Fock matrix (F + Delta')
    for(auto iVir = this->nOccB_; iVir < this->nBasis_; iVir++){
      (*this->fockB_)(iVir, iVir) += deltaBPrime(iVir, iVir);	//Sums the diagonal elements of the vir-vir blocks only
    }
  }
}
*/

template<typename T>
void SingleSlater<T>::levelShift2(){
  if(this->levelShiftParam <= 0.0) return;
/*
  cout << "LEVEL Shifting" << endl;
  T* FockA, *FockB;
  T* MOAVir, *MOBVir;
  T MOAMu, MOBMu;
//T* FockA = this->fockA_->data();
  FockA = this->fockOrthoA_->data();
  if(this->nTCS_ == 1 && !this->isClosedShell) 
    FockB = this->fockOrthoB_->data();
  
  for(auto iVir = this->nOccA_; iVir < this->nBasis_; iVir++){
    MOAVir = this->moA_->data() + iVir * this->nBasis_;
    if(this->nTCS_ == 1 && !this->isClosedShell) 
      MOBVir = this->moB_->data() + iVir * this->nBasis_;

    for(auto mu = 0; mu < this->nBasis_; mu++){

      MOAMu = this->levelShiftParam_ * MOAVir[mu];
      for(auto nu = 0; nu < this->nBasis_; nu++)
        FockA[nu + mu*this->nBasis_] += MOAMu * MOAVir[nu];

      if(this->nTCS_ == 1 && !this->isClosedShell){ 
        MOBMu = this->levelShiftParam_ * MOBVir[mu];
        for(auto nu = 0; nu < this->nBasis_; nu++)
          FockB[nu + mu*this->nBasis_] += MOBMu * MOBVir[nu];
      }
    
    }
  }
*/
  double HOMO_LUMO_GAP = 
    (*this->epsA_)(this->nO_) - (*this->epsA_)(this->nO_-1);

  if( HOMO_LUMO_GAP < 1e-6 ) {
    this->levelShiftParam = 0.1;
  } else {
    this->levelShiftParam = HOMO_LUMO_GAP * 2;
  }

  T *MOAVir, *MOBVir;
  T MOAMu, MOBMu;

  T *FockScalar, *FockMz, *FockMy, *FockMx;
  FockScalar = this->fockOrthoScalar_->data();
  if(this->nTCS_ == 2 or !this->isClosedShell)
    FockMz = this->fockOrthoMz_->data();
  if(this->nTCS_ == 2){
    FockMy = this->fockOrthoMy_->data();
    FockMx = this->fockOrthoMx_->data();
  }

  int NO = this->nTCS_ == 2 ? this->nO_ : this->nOA_;
  for(auto iVir = NO; iVir < this->nTCS_ * this->nBasis_; iVir++){
    MOAVir = this->moA_->data() + iVir * this->nTCS_ * this->nBasis_;
    if(this->nTCS_ == 1 && this->isClosedShell) 
      MOBVir = this->moA_->data() + iVir * this->nBasis_;
    else if(this->nTCS_ == 1)
      MOBVir = this->moB_->data() + iVir * this->nBasis_;
    else
      MOBVir = this->moA_->data() + 1 + iVir * this->nTCS_ * this->nBasis_;

    for(auto mu = 0; mu < this->nTCS_ * this->nBasis_; mu += this->nTCS_){
      MOAMu = this->levelShiftParam * MOAVir[mu];
      MOBMu = this->levelShiftParam * MOBVir[mu];

      // Scalar
      for(auto nu = 0; nu < this->nTCS_ * this->nBasis_; nu += this->nTCS_){
        FockScalar[nu/this->nTCS_ + (mu/this->nTCS_) * this->nBasis_] +=
          MOAMu * std::conj(MOAVir[nu]) + MOBMu * std::conj(MOBVir[nu]);
      }
      if(this->nTCS_ == 2 or !this->isClosedShell){
        // Mz 
        for(auto nu = 0; nu < this->nTCS_ * this->nBasis_; nu += this->nTCS_){
          FockMz[nu/this->nTCS_ + (mu/this->nTCS_) * this->nBasis_] +=
            MOAMu * std::conj(MOAVir[nu]) - MOBMu * std::conj(MOBVir[nu]);
        }
      }
      if(this->nTCS_ == 2 ){
        // Mx / My
        for(auto nu = 0; nu < this->nTCS_ * this->nBasis_; nu += this->nTCS_){
          FockMx[nu/this->nTCS_ + (mu/this->nTCS_) * this->nBasis_] +=
            MOAMu * std::conj(MOBVir[nu]) + MOBMu * std::conj(MOAVir[nu]);
          FockMy[nu/this->nTCS_ + (mu/this->nTCS_) * this->nBasis_] +=
            ComplexScale<T>() * (
            MOAMu * std::conj(MOBVir[nu]) - MOBMu * std::conj(MOAVir[nu])
            );
        }
      }
    }
    
  }
  
}
// JJRADLER
//
//
//

