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
/***********************
 * Form Density Matrix *
 ***********************/
template<typename T>
void SingleSlater<T>::formDensity(){
  if(getRank() == 0) {
    if(!this->haveMO)
      CErr("No MO coefficients available to form one-particle density matrix!",
           this->fileio_->out);
 
    if(this->nTCS_ == 2){
      auto nOcc = this->nOccA_ + this->nOccB_;
      (*this->onePDMA_) = 
        this->moA_->block(0,0,this->nTCS_*this->nBasis_,nOcc)*
        this->moA_->block(0,0,this->nTCS_*this->nBasis_,nOcc).adjoint();
    } else {
      (*this->onePDMA_) = 
        this->moA_->block(0,0,this->nBasis_,this->nOccA_)*
        this->moA_->block(0,0,this->nBasis_,this->nOccA_).adjoint();
      // D(a) is actually total D for RHF
      if(this->isClosedShell) (*this->onePDMA_) *= math.two;
      else {
        (*this->onePDMB_) = 
          this->moB_->block(0,0,this->nBasis_,this->nOccB_)*
          this->moB_->block(0,0,this->nBasis_,this->nOccB_).adjoint();
      }
    }
    if(this->printLevel_ >= 2) this->printDensity();
  }
  this->haveDensity = true;
  this->mpiBCastDensity();
}

template<typename T>
void SingleSlater<T>::formFP(){
  // FP(S) = F(S)P(S)
  if(this->nTCS_ == 1 and this->isClosedShell) {
    (*this->NBSqScratch_) = (*this->fockOrthoA_) * (*this->onePDMOrthoA_);
  } else {
    (*this->NBSqScratch_) = 
      (*this->fockOrthoScalar_) * (*this->onePDMOrthoScalar_);
  } 

  // FP(S) += F(z)P(z)
  if(this->nTCS_ == 2 or !this->isClosedShell)
    (*this->NBSqScratch_) += (*this->fockOrthoMz_) * (*this->onePDMOrthoMz_);

  // FP(S) += F(y)P(y) + F(x)P(x)
  if(this->nTCS_ == 2) {
    (*this->NBSqScratch_) += (*this->fockOrthoMy_) * (*this->onePDMOrthoMy_);
    (*this->NBSqScratch_) += (*this->fockOrthoMx_) * (*this->onePDMOrthoMx_);
  }

  this->FPScalar_->write(this->NBSqScratch_->data(),H5PredType<T>());

  
  if(this->nTCS_ == 2 or !this->isClosedShell) {
    // FP(z) = F(S)P(z) + F(z)P(S)
    (*this->NBSqScratch_) = 
      (*this->fockOrthoScalar_) * (*this->onePDMOrthoMz_);
    (*this->NBSqScratch_) += 
      (*this->fockOrthoMz_) * (*this->onePDMOrthoScalar_);

    // FP(z) += i * (F(x)P(y) - F(y)P(x))
    if(this->nTCS_ == 2) {
      (*this->NBSqScratch_) += DIISComplexScale<T>() * 
        (*this->fockOrthoMx_) * (*this->onePDMOrthoMy_);
      (*this->NBSqScratch_) -= DIISComplexScale<T>() * 
        (*this->fockOrthoMy_) * (*this->onePDMOrthoMx_);
    }
    this->FPMz_->write(this->NBSqScratch_->data(),H5PredType<T>());
  }


  if(this->nTCS_ == 2) {
    // FP(y) = F(S)P(y) + F(y)P(S)
    (*this->NBSqScratch_) = 
      (*this->fockOrthoScalar_) * (*this->onePDMOrthoMy_);
    (*this->NBSqScratch_) += 
      (*this->fockOrthoMy_) * (*this->onePDMOrthoScalar_);

    // FP(y) += i * (F(z)P(x) - F(x)P(z))
    (*this->NBSqScratch_) += DIISComplexScale<T>() * 
      (*this->fockOrthoMz_) * (*this->onePDMOrthoMx_);
    (*this->NBSqScratch_) -= DIISComplexScale<T>() * 
      (*this->fockOrthoMx_) * (*this->onePDMOrthoMz_);

    this->FPMy_->write(this->NBSqScratch_->data(),H5PredType<T>());

    // FP(x) = F(S)P(x) + F(x)P(S)
    (*this->NBSqScratch_) = 
      (*this->fockOrthoScalar_) * (*this->onePDMOrthoMx_);
    (*this->NBSqScratch_) += 
      (*this->fockOrthoMx_) * (*this->onePDMOrthoScalar_);

    // FP(x) += i * (F(y)P(z) - F(z)P(y))
    (*this->NBSqScratch_) += DIISComplexScale<T>() * 
      (*this->fockOrthoMy_) * (*this->onePDMOrthoMz_);
    (*this->NBSqScratch_) -= DIISComplexScale<T>() * 
      (*this->fockOrthoMz_) * (*this->onePDMOrthoMy_);

    this->FPMx_->write(this->NBSqScratch_->data(),H5PredType<T>());
  }
};
