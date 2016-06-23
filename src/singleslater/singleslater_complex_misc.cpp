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
namespace ChronusQ {
template<>
template<>
SingleSlater<dcomplex>::SingleSlater(SingleSlater<double> * other) :
  Quantum<dcomplex>::Quantum<dcomplex>(dynamic_cast<Quantum<double>&>(*other)){
     

    this->nBasis_ = other->nBasis();
    this->nTT_    = other->nTT();
    this->nAE_    = other->nAE();
    this->nBE_    = other->nBE(); 
    this->nOccA_  = other->nOccA();
    this->nOccB_  = other->nOccB();
    this->nVirA_  = other->nVirA();
    this->nVirB_  = other->nVirB();
    this->multip_   = other->multip();
    this->energyNuclei = other->energyNuclei;
    this->Ref_    = other->Ref();
    this->haveDensity = true;
    this->haveMO	    = true;
    this->havePT      = true;
    this->printLevel_ = other->printLevel();
    this->doDIIS = other->doDIIS;
    this->isHF   = other->isHF;
    this->isDFT  = other->isDFT;
    this->guess_ = other->guess();
    this->elecField_   = (other->elecField());
    this->basisset_    = other->basisset();    
    this->molecule_    = other->molecule();
    this->fileio_      = other->fileio();
    this->aointegrals_ = other->aointegrals();

    auto NB = this->nBasis_*this->nTCS_;
    auto NBSq = NB*NB;
    this->allocOp();

    this->fockA_->real()       = *other->fockA();
    this->moA_->real()         = *other->moA();
    this->PTA_->real()         = *other->PTA();

    if(!this->isClosedShell || this->nTCS_ == 2){
      this->fockB_->real()       = *other->fockB();
      this->moB_->real()         = *other->moB();
      this->PTB_->real()         = *other->PTB();
    }

}

template<>
void SingleSlater<dcomplex>::getAlgebraicField(){ 
  this->algebraicField_      = "Complex";
  this->algebraicFieldShort_ = "\u2102";
}

template<>
void SingleSlater<dcomplex>::writeSCFFiles(){

  this->fileio_->alphaSCFDen->write(this->onePDMA_->data(),
    *(this->fileio_->complexType)
  );
  this->fileio_->alphaMO->write(this->moA_->data(),
    *(this->fileio_->complexType)
  );
  if(!this->isClosedShell && this->Ref_ != TCS){
    this->fileio_->betaSCFDen->write(this->onePDMB_->data(),
      *(this->fileio_->complexType)
    );
    this->fileio_->betaMO->write(this->moB_->data(),
      *(this->fileio_->complexType)
    );
  }
}

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
