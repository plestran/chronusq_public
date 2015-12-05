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

template<typename T>
void QuasiNewton2<T>::allocScr(){

  // Linear Dimensions
  auto N      = this->qnObj_->nSingleDim();
  auto NMSS   = N*this->maxSubSpace_;
  auto MSSMSS = this->maxSubSpace_*this->maxSubSpace_;

  this->TRMem_       = new T[NMSS];
  this->SigmaRMem_   = new T[NMSS]; 
  this->XTSigmaRMem_ = new T[MSSMSS];
  this->ResRMem_     = new T[NMSS];
  this->URMem_       = new T[NMSS];

  if(this->matrixType_ == HERMETIAN_GEP){
    this->RhoRMem_   = new T[NMSS];
    this->XTRhoRMem_ = new T[MSSMSS];
  }

  if(this->qnObj_->needsLeft()){
    this->TLMem_       = new T[NMSS];
    this->SigmaLMem_   = new T[NMSS]; 
    this->XTSigmaLMem_ = new T[MSSMSS];
    this->ResLMem_     = new T[NMSS];
    this->ULMem_       = new T[NMSS];
 
    if(this->matrixType_ == HERMETIAN_GEP){
      this->RhoLMem_   = new T[NMSS];
      this->XTRhoLMem_ = new T[MSSMSS];
    }
  }

  if(this->problemType_ == DIAGONALIZATION){
    if(this->matrixType_ == HERMETIAN)
      this->LWORK = 3 * N;
    else
      this->LWORK = 4 * N;


    this->WORK   = new T[this->LWORK];
  }

  this->allocScrSpecial();
};
