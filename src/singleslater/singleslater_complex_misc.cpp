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
/*
template<>
template<>
SingleSlater<dcomplex>::SingleSlater(SingleSlater<double> * other) :
  WaveFunction<dcomplex>::WaveFunction<dcomplex>(
    dynamic_cast<WaveFunction<double>&>(*other)),
  Ref_    ( other->Ref() ),
  printLevel_  ( other->printLevel() ),
  doDIIS  ( other->doDIIS ),
  isHF    ( other->isHF ),
  isDFT   ( other->isDFT ),
  dftFunctionals_( other->dftFunctionals_ ),
  weightScheme_( other->weightScheme()),
  dftGrid_( other->dftGrid()),
  nRadDFTGridPts_( other->nRadGridPts() ),
  nAngDFTGridPts_( other->nAngGridPts() ),
  xHF_( other->xHF() ),
  isGGA(other->isGGA),
  screenVxc( other->screenVxc),
  epsScreen( other->epsScreen),
  epsConv( other->epsConv),
  maxiter( other->maxiter),
  guess_  ( other->guess() ),
  elecField_   ( other->elecField() ) {

    this->alloc();

    for(auto F : this->dftFunctionals_) cout << F->name << endl;
    cout << this->xHF_ << endl;
    for(auto iF = 0; iF < this->fock_.size(); iF++){
      this->fock_[iF]->real() = *other->fock()[iF];
      this->PT_[iF]->real()   = *other->PT()[iF];

      // Copy over the orthonormal density for RT calculations
      this->onePDMOrtho_[iF]->real() = *other->onePDMOrtho()[iF];
    }


}
*/

template<>
void SingleSlater<dcomplex>::getAlgebraicField(){ 
  this->algebraicField_      = "Complex";
  this->algebraicFieldShort_ = "\u2102";
}

/*
template<>
void SingleSlater<dcomplex>::writeSCFFiles(){

  this->fileio_->alphaSCFDen->write(this->onePDMA_->data(),
    *(this->fileio_->complexType)
  );
  this->fileio_->alphaMO->write(this->moA_->data(),
    *(this->fileio_->complexType)
  );
  if(!this->isClosedShell && this->nTCS_ == 1){
    this->fileio_->betaSCFDen->write(this->onePDMB_->data(),
      *(this->fileio_->complexType)
    );
    this->fileio_->betaMO->write(this->moB_->data(),
      *(this->fileio_->complexType)
    );
  }
}
*/

template<>
void SingleSlater<dcomplex>::fixPhase(){
  // FIXME: Do nothing for now
};

template<>
void SingleSlater<dcomplex>::backTransformMOs(){
  if(this->nTCS_ == 1) {
    this->NBSqScratch_->real() = 
      (*this->aointegrals_->ortho1_) * this->moA_->real();
    this->NBSqScratch_->imag() = 
      (*this->aointegrals_->ortho1_) * this->moA_->imag();

    (*this->moA_) = (*this->NBSqScratch_);

    if(!this->isClosedShell){
      this->NBSqScratch_->real() = 
        (*this->aointegrals_->ortho1_) * this->moB_->real();
      this->NBSqScratch_->imag() = 
        (*this->aointegrals_->ortho1_) * this->moB_->imag();
      (*this->moB_) = (*this->NBSqScratch_);
    }
  } else {
    Eigen::Map<ComplexMatrix,0,Eigen::Stride<Dynamic,Dynamic> >
      MOA(this->moA_->data(),this->nBasis_,this->nTCS_*this->nBasis_,
          Eigen::Stride<Dynamic,Dynamic>(this->nTCS_*this->nBasis_,2));
    Eigen::Map<ComplexMatrix,0,Eigen::Stride<Dynamic,Dynamic> >
      MOB(this->moA_->data()+1,this->nBasis_,this->nTCS_*this->nBasis_,
          Eigen::Stride<Dynamic,Dynamic>(this->nTCS_*this->nBasis_,2));

    ComplexMap SCRATCH1(this->memManager_->malloc<dcomplex>(this->nBasis_*
          this->nBasis_*this->nTCS_),this->nBasis_,this->nTCS_*this->nBasis_);
    ComplexMap SCRATCH2(this->memManager_->malloc<dcomplex>(this->nBasis_*
          this->nBasis_*this->nTCS_),this->nBasis_,this->nTCS_*this->nBasis_);

    SCRATCH1 = MOA;
    SCRATCH2.real() =
      (*this->aointegrals_->ortho1_) * SCRATCH1.real();
    SCRATCH2.imag() =
      (*this->aointegrals_->ortho1_) * SCRATCH1.imag();
    MOA = SCRATCH2;

    SCRATCH1 = MOB;
    SCRATCH2.real() =
      (*this->aointegrals_->ortho1_) * SCRATCH1.real();
    SCRATCH2.imag() =
      (*this->aointegrals_->ortho1_) * SCRATCH1.imag();
    MOB = SCRATCH2;

    this->memManager_->free(SCRATCH1.data(),
        this->nBasis_*this->nBasis_*this->nTCS_);
    this->memManager_->free(SCRATCH2.data(),
        this->nBasis_*this->nBasis_*this->nTCS_);
    
  }
};
} // Namespace ChronusQ
