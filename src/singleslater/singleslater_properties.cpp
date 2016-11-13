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
#include <singleslater.h>
namespace ChronusQ {
template<>
void SingleSlater<double>::computeSSq(){
  (*this->NBSqScratch_)  = (*this->onePDMMz_) * (*this->aointegrals_->overlap_);
  (*this->NBSqScratch2_) = 
    (*this->aointegrals_->overlap_) * (*this->NBSqScratch_);

  double M = this->computeProperty<double,MZ>(*this->aointegrals_->overlap_);
  this->Ssq_ = M * M + 2 * 
    this->computeProperty<double,MZ>(*this->NBSqScratch2_);

  if(this->nTCS_ == 2) {
    (*this->NBSqScratch_)  = 
      (*this->onePDMMx_) * (*this->aointegrals_->overlap_);
    (*this->NBSqScratch2_) = 
      (*this->aointegrals_->overlap_) * (*this->NBSqScratch_);

    M = this->computeProperty<double,MX>(*this->aointegrals_->overlap_);
    this->Ssq_ += M * M + 2 * 
      this->computeProperty<double,MX>(*this->NBSqScratch2_);
  }
  this->Ssq_ *= 0.25;
};
template<>
void SingleSlater<dcomplex>::computeSSq(){
  this->NBSqScratch_->real() = 
    this->onePDMMz_->real() * (*this->aointegrals_->overlap_);
  this->NBSqScratch_->imag() = 
    this->onePDMMz_->imag() * (*this->aointegrals_->overlap_);

  this->NBSqScratch2_->real() = 
    (*this->aointegrals_->overlap_) * this->NBSqScratch_->real();
  this->NBSqScratch2_->imag() = 
    (*this->aointegrals_->overlap_) * this->NBSqScratch_->imag();

  double M = this->computeProperty<double,MZ>(*this->aointegrals_->overlap_);
  this->Ssq_ = M * M + 2 * 
    this->computeProperty<double,MZ>(*this->NBSqScratch2_);

  if(this->nTCS_ == 2) {
    this->NBSqScratch_->real() = 
      this->onePDMMx_->real() * (*this->aointegrals_->overlap_);
    this->NBSqScratch_->imag() = 
      this->onePDMMx_->imag() * (*this->aointegrals_->overlap_);
 
    this->NBSqScratch2_->real() = 
      (*this->aointegrals_->overlap_) * this->NBSqScratch_->real();
    this->NBSqScratch2_->imag() = 
      (*this->aointegrals_->overlap_) * this->NBSqScratch_->imag();

    M = this->computeProperty<double,MX>(*this->aointegrals_->overlap_);
    this->Ssq_ += M * M + 2 * 
      this->computeProperty<double,MX>(*this->NBSqScratch2_);

    this->NBSqScratch_->real() = 
      this->onePDMMy_->real() * (*this->aointegrals_->overlap_);
    this->NBSqScratch_->imag() = 
      this->onePDMMy_->imag() * (*this->aointegrals_->overlap_);
 
    this->NBSqScratch2_->real() = 
      (*this->aointegrals_->overlap_) * this->NBSqScratch_->real();
    this->NBSqScratch2_->imag() = 
      (*this->aointegrals_->overlap_) * this->NBSqScratch_->imag();

    M = this->computeProperty<double,MY>(*this->aointegrals_->overlap_);
    this->Ssq_ += M * M + 2 * 
      this->computeProperty<double,MY>(*this->NBSqScratch2_);
  }
  this->Ssq_ *= 0.25;
};
};
