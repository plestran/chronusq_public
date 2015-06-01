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
    this->RHF_    = other->RHF_;
    this->nOccA_  = other->nOccA_;
    this->nOccB_  = other->nOccB_;
    this->nVirA_  = other->nVirA_;
    this->nVirB_  = other->nVirB_;
    this->spin_   = other->spin_;
    // Hardcoded for Libint route
    this->densityA_           = std::unique_ptr<ComplexMatrix>(new ComplexMatrix(*other->densityA_));
    this->fockA_              = std::unique_ptr<ComplexMatrix>(new ComplexMatrix(*other->fockA_));
    this->moA_                = std::unique_ptr<ComplexMatrix>(new ComplexMatrix(*other->moA_));
    this->PTA_                = std::unique_ptr<ComplexMatrix>(new ComplexMatrix(*other->PTA_));
    if(!this->RHF_){
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
    this->RHF_    = other->RHF();
    this->nOccA_  = other->nOccA();
    this->nOccB_  = other->nOccB();
    this->nVirA_  = other->nVirA();
    this->nVirB_  = other->nVirB();
    this->spin_   = other->spin();
    // Hardcoded for Libint route
    this->densityA_           = std::unique_ptr<ComplexMatrix>(new ComplexMatrix(this->nBasis_,this->nBasis_));
    this->fockA_              = std::unique_ptr<ComplexMatrix>(new ComplexMatrix(this->nBasis_,this->nBasis_));
    this->moA_                = std::unique_ptr<ComplexMatrix>(new ComplexMatrix(this->nBasis_,this->nBasis_));
    this->PTA_                = std::unique_ptr<ComplexMatrix>(new ComplexMatrix(this->nBasis_,this->nBasis_));
    this->densityA_->real()    = *other->densityA();
    this->fockA_->real()       = *other->fockA();
    this->moA_->real()         = *other->moA();
    this->PTA_->real()         = *other->PTA();
    if(!this->RHF_){
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
} // Namespace ChronusQ
