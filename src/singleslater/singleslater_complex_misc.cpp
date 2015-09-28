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
namespace ChronusQ {
template<>
template<>
SingleSlater<dcomplex>::SingleSlater(SingleSlater<dcomplex> * other){
    this->nBasis_ = other->nBasis_;
    this->nTT_    = other->nTT_;
    this->nAE_    = other->nAE_;
    this->nBE_    = other->nBE_; 
    this->Ref_    = other->Ref();
    this->nOccA_  = other->nOccA_;
    this->nOccB_  = other->nOccB_;
    this->nVirA_  = other->nVirA_;
    this->nVirB_  = other->nVirB_;
    this->multip_   = other->multip_;
    this->energyNuclei = other->energyNuclei;
    this->haveDensity = true;
    this->haveMO	    = true;
    this->havePT      = true;
    this->isClosedShell = other->isClosedShell;
    // Hardcoded for Libint route
    this->densityA_           = std::unique_ptr<ComplexMatrix>(new ComplexMatrix(*other->densityA_));
    this->fockA_              = std::unique_ptr<ComplexMatrix>(new ComplexMatrix(*other->fockA_));
    this->moA_                = std::unique_ptr<ComplexMatrix>(new ComplexMatrix(*other->moA_));
    this->PTA_                = std::unique_ptr<ComplexMatrix>(new ComplexMatrix(*other->PTA_));
    if(this->Ref_ != RHF){
      this->densityB_           = std::unique_ptr<ComplexMatrix>(new ComplexMatrix(*other->densityB_));
      this->fockB_              = std::unique_ptr<ComplexMatrix>(new ComplexMatrix(*other->fockB_));
      this->moB_                = std::unique_ptr<ComplexMatrix>(new ComplexMatrix(*other->moB_));
      this->PTB_                = std::unique_ptr<ComplexMatrix>(new ComplexMatrix(*other->PTB_));
    }
    this->dipole_             = std::unique_ptr<RealMatrix>(new RealMatrix(*other->dipole_));
    this->quadpole_           = std::unique_ptr<RealMatrix>(new RealMatrix(*other->quadpole_));
    this->tracelessQuadpole_  = std::unique_ptr<RealMatrix>(new RealMatrix(*other->tracelessQuadpole_));
    this->octpole_            = std::unique_ptr<RealTensor3d>(new RealTensor3d(*other->octpole_));
    this->basisset_    = other->basisset_;    
    this->molecule_    = other->molecule_;
    this->fileio_      = other->fileio_;
    this->controls_    = other->controls_;
    this->aointegrals_ = other->aointegrals_;
}
template<>
template<>
SingleSlater<dcomplex>::SingleSlater(SingleSlater<double> * other){
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
    this->isClosedShell = other->isClosedShell;
    // Hardcoded for Libint route
    this->densityA_           = std::unique_ptr<ComplexMatrix>(new ComplexMatrix(this->nBasis_,this->nBasis_));
    this->fockA_              = std::unique_ptr<ComplexMatrix>(new ComplexMatrix(this->nBasis_,this->nBasis_));
    this->moA_                = std::unique_ptr<ComplexMatrix>(new ComplexMatrix(this->nBasis_,this->nBasis_));
    this->PTA_                = std::unique_ptr<ComplexMatrix>(new ComplexMatrix(this->nBasis_,this->nBasis_));
    this->densityA_->real()    = *other->densityA();
    this->fockA_->real()       = *other->fockA();
    this->moA_->real()         = *other->moA();
    this->PTA_->real()         = *other->PTA();
    if(this->Ref_ != RHF){
      this->densityB_           = std::unique_ptr<ComplexMatrix>(new ComplexMatrix(this->nBasis_,this->nBasis_));
      this->fockB_              = std::unique_ptr<ComplexMatrix>(new ComplexMatrix(this->nBasis_,this->nBasis_));
      this->moB_                = std::unique_ptr<ComplexMatrix>(new ComplexMatrix(this->nBasis_,this->nBasis_));
      this->PTB_                = std::unique_ptr<ComplexMatrix>(new ComplexMatrix(this->nBasis_,this->nBasis_));
      this->densityB_->real()    = *other->densityB();
      this->fockB_->real()       = *other->fockB();
      this->moB_->real()         = *other->moB();
      this->PTB_->real()         = *other->PTB();
    }
    this->dipole_             = std::unique_ptr<RealMatrix>(new RealMatrix(*other->dipole()));
    this->quadpole_           = std::unique_ptr<RealMatrix>(new RealMatrix(*other->quadpole()));
    this->tracelessQuadpole_  = std::unique_ptr<RealMatrix>(new RealMatrix(*other->tracelessQuadpole()));
    this->octpole_            = std::unique_ptr<RealTensor3d>(new RealTensor3d(*other->octpole()));
    this->basisset_    = other->basisset();    
    this->molecule_    = other->molecule();
    this->fileio_      = other->fileio();
    this->controls_    = other->controls();
    this->aointegrals_ = other->aointegrals();
}
/************************
 * Compute Total Energy *
 ************************/
template<>
void SingleSlater<dcomplex>::computeEnergy(){
/*
  if(this->Ref_ != TCS)
    this->energyOneE = (*this->aointegrals_->oneE_).frobInner(this->densityA_->conjugate());
  else {
    this->energyOneE = 0.0;
    for(auto I = 0, i = 0; i < this->nBasis_; I += 2, i++)    
    for(auto J = 0, j = 0; j < this->nBasis_; J += 2, j++){
      this->energyOneE += 
        this->densityA_->conjugate()(I,J)*(*this->aointegrals_->oneE_)(i,j) + 
        this->densityA_->conjugate()(I+1,J+1)*(*this->aointegrals_->oneE_)(i,j);
    } 
  }
*/
  this->energyOneE = (*this->aointegrals_->oneE_).frobInner(this->densityA_->real());
  this->energyTwoE = 0.5*(*this->PTA_).frobInner(this->densityA_->conjugate()).real();
  if(!this->isClosedShell && this->Ref_ != TCS){
    this->energyOneE += (*this->aointegrals_->oneE_).frobInner(this->densityB_->real());
    this->energyTwoE += 0.5*(*this->PTB_).frobInner(this->densityB_->conjugate()).real();
  }
  this->totalEnergy= this->energyOneE + this->energyTwoE + this->energyNuclei;
  this->printEnergy();
};
} // Namespace ChronusQ
