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
/*
  if(getRank() == 0) {
 
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
*/
  if(this->nTCS_ == 1) {
    // Store Pa in Ps
    this->onePDMScalar_->noalias() = 
      this->moA_->block(0,0,this->nBasis_,this->nOccA_)*
      this->moA_->block(0,0,this->nBasis_,this->nOccA_).adjoint();
    if(!this->isClosedShell) {
      // Store Pb in Scratch
      this->NBSqScratch_->noalias() = 
        this->moB_->block(0,0,this->nBasis_,this->nOccB_)*
        this->moB_->block(0,0,this->nBasis_,this->nOccB_).adjoint();
      // Overwrite Pz with Pa - Pb
      (*this->onePDMMz_) =     (*this->onePDMScalar_) - (*this->NBSqScratch_);
      // Overwrite Ps with Pa + Pb
      (*this->onePDMScalar_) = (*this->onePDMScalar_) + (*this->NBSqScratch_);
    } else {
      // Factor of 2 for scalar
      (*this->onePDMScalar_) *= 2;
    }
  } else {
    auto NE = this->nOccA_ + this->nOccB_;
    this->NBTSqScratch_->noalias() = 
      this->moA_->block(0,0,this->nTCS_*this->nBasis_,NE)*
      this->moA_->block(0,0,this->nTCS_*this->nBasis_,NE).adjoint();
    std::vector<std::reference_wrapper<TMap>> scatter;
    for(auto iD = 0; iD < this->onePDM_.size(); iD++)
      scatter.emplace_back(*this->onePDM_[iD]);
    this->spinScatter((*this->NBTSqScratch_),scatter); 
  }
}

template<typename T>
void SingleSlater<T>::formFP(){
  // FP(S) = F(S)P(S)
  (*this->NBSqScratch_) = 
    (*this->fockOrthoScalar_) * (*this->onePDMOrthoScalar_);

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

template<typename T>
void SingleSlater<T>::fockCUHF() {
  TMap P(this->PNOMem_,this->nBasis_,this->nBasis_);
  TMap DelF(this->delFMem_,this->nBasis_,this->nBasis_);
  TMap Lambda(this->lambdaMem_,this->nBasis_,this->nBasis_);

  int activeSpace  = this->molecule_->multip() - 1;
  int coreSpace    = (this->molecule_->nTotalE() - activeSpace) / 2;
  int virtualSpace = this->nBasis_ - coreSpace - activeSpace;

  // DelF = X * (F(A) - F(B)) * X

  this->aointegrals_->Ortho1Trans(*this->fockMz_,DelF);
  (*this->fockMz_) *= 0.5;

  // DelF = C(NO)^\dagger * DelF * C(NO) (Natural Orbitals)
  (*this->NBSqScratch_) = P.transpose() * DelF;
  DelF = (*this->NBSqScratch_) * P;

  Lambda.setZero();
  for(auto i = activeSpace + coreSpace; i < this->nBasis_; i++)
  for(auto j = 0                      ; j < coreSpace    ; j++){
    Lambda(i,j) = -DelF(i,j);
    Lambda(j,i) = -DelF(j,i);
  }

  (*this->NBSqScratch_) = P * Lambda;
  Lambda = (*this->NBSqScratch_) * P.transpose();

  this->aointegrals_->Ortho2Trans(Lambda,Lambda);

  (*this->fockMz_) += 2*Lambda;
};

