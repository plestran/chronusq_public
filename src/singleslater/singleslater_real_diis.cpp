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
template<>
void SingleSlater<double>::GenDComm(int iter){
  RealMap ErrA(
    this->ErrorAlphaMem_ + (iter % (this->nDIISExtrap_-1)) * this->lenF_,
    this->nTCS_*this->nBasis_,this->nTCS_*this->nBasis_
  );

  ErrA = (*this->fockA_) * (*this->onePDMA_) * (*this->aointegrals_->overlap_);
  ErrA -= (*this->aointegrals_->overlap_) * (*this->onePDMA_) * (*this->fockA_);
  if(!this->isClosedShell && this->nTCS_ == 1){
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
