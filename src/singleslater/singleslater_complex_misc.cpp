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
SingleSlater<dcomplex>::SingleSlater(SingleSlater<dcomplex> * other){
    this->nBasis_ = other->nBasis_;
    this->nTCS_   = other->nTCS_;
    this->nTT_    = other->nTT_;
    this->nShell_ = other->nShell_;
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
    this->printLevel_ = other->printLevel_;
    this->maxMultipole_ = other->maxMultipole_;
    this->doDIIS = other->doDIIS;
    this->isHF   = other->isHF;
    this->isDFT  = other->isDFT;
    this->guess_ = other->guess_;

    auto NTCSxNBASIS = this->nBasis_*this->nTCS_;

    // Hardcoded for Libint route
    this->densityA_ = std::unique_ptr<ComplexMatrix>(
      new ComplexMatrix(*other->densityA_)
    );
    if(getRank() == 0) {
      this->fockA_    = std::unique_ptr<ComplexMatrix>(
        new ComplexMatrix(*other->fockA_)
      );
      this->moA_      = std::unique_ptr<ComplexMatrix>(
        new ComplexMatrix(*other->moA_)
      );
      this->PTA_      = std::unique_ptr<ComplexMatrix>(
        new ComplexMatrix(*other->PTA_)
      );
    }
    if(this->Ref_ != isClosedShell && this->Ref_ != TCS ){
      this->densityB_ = std::unique_ptr<ComplexMatrix>(
        new ComplexMatrix(*other->densityB_)
      );
      if(getRank() == 0) {
        this->fockB_    = std::unique_ptr<ComplexMatrix>(
          new ComplexMatrix(*other->fockB_)
        );
        this->moB_      = std::unique_ptr<ComplexMatrix>(
          new ComplexMatrix(*other->moB_)
        );
        this->PTB_      = std::unique_ptr<ComplexMatrix>(
          new ComplexMatrix(*other->PTB_)
        );
      }
    }
    this->dipole_             = std::unique_ptr<RealMatrix>(
      new RealMatrix(*other->dipole_)
    );
    this->quadpole_           = std::unique_ptr<RealMatrix>(
      new RealMatrix(*other->quadpole_)
    );
    this->tracelessQuadpole_  = std::unique_ptr<RealMatrix>(
      new RealMatrix(*other->tracelessQuadpole_)
    );
    this->octpole_            = std::unique_ptr<RealTensor3d>(
      new RealTensor3d(*other->octpole_)
    );

    this->elecField_   = other->elecField_;
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
    this->nTCS_   = other->nTCS();
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
    this->printLevel_ = other->printLevel();
    this->maxMultipole_ = other->maxMultipole();
    this->doDIIS = other->doDIIS;
    this->isHF   = other->isHF;
    this->isDFT  = other->isDFT;
    this->guess_ = other->guess();

    auto NTCSxNBASIS = this->nBasis_*this->nTCS_;

    this->densityA_ = std::unique_ptr<ComplexMatrix>(
      new ComplexMatrix(NTCSxNBASIS,NTCSxNBASIS)
    );

    if(getRank() == 0) {
      this->fockA_    = std::unique_ptr<ComplexMatrix>(
        new ComplexMatrix(NTCSxNBASIS,NTCSxNBASIS)
      );
      this->moA_      = std::unique_ptr<ComplexMatrix>(
        new ComplexMatrix(NTCSxNBASIS,NTCSxNBASIS)
      );
      this->PTA_      = std::unique_ptr<ComplexMatrix>(
        new ComplexMatrix(NTCSxNBASIS,NTCSxNBASIS)
      );
    }

    this->densityA_->real()    = *other->densityA();
    if(getRank() == 0) {
      this->fockA_->real()       = *other->fockA();
      this->moA_->real()         = *other->moA();
      this->PTA_->real()         = *other->PTA();
    }
    if(this->Ref_ != isClosedShell && this->Ref_ != TCS ){
      this->densityB_           = std::unique_ptr<ComplexMatrix>(
        new ComplexMatrix(NTCSxNBASIS,NTCSxNBASIS)
      );
     
      if(getRank() == 0) {
        this->fockB_  = std::unique_ptr<ComplexMatrix>(
          new ComplexMatrix(NTCSxNBASIS,NTCSxNBASIS)
        );
        this->moB_    = std::unique_ptr<ComplexMatrix>(
          new ComplexMatrix(NTCSxNBASIS,NTCSxNBASIS)
        );
        this->PTB_    = std::unique_ptr<ComplexMatrix>(
          new ComplexMatrix(NTCSxNBASIS,NTCSxNBASIS)
        );
      }

      this->densityB_->real()    = *other->densityB();
      if(getRank() == 0) {
        this->fockB_->real()       = *other->fockB();
        this->moB_->real()         = *other->moB();
        this->PTB_->real()         = *other->PTB();
      }
    }
    this->dipole_ = std::unique_ptr<RealMatrix>(
      new RealMatrix(*other->dipole())
    );
    this->quadpole_ = std::unique_ptr<RealMatrix>(
      new RealMatrix(*other->quadpole())
    );
    this->tracelessQuadpole_  = std::unique_ptr<RealMatrix>(
      new RealMatrix(*other->tracelessQuadpole())
    );
    this->octpole_ = std::unique_ptr<RealTensor3d>(
      new RealTensor3d(*other->octpole())
     );

    this->elecField_   = (other->elecField());
    this->basisset_    = other->basisset();    
    this->molecule_    = other->molecule();
    this->fileio_      = other->fileio();
    this->controls_    = other->controls();
    this->aointegrals_ = other->aointegrals();
}
/************************
 * Compute Total Energy *
 ************************/
/*
template<>
void SingleSlater<dcomplex>::computeEnergy(){
  if(getRank() == 0) {
    this->energyOneE = 
      (*this->aointegrals_->oneE_).frobInner(this->densityA_->real());
    this->energyTwoE = 
      0.5*(*this->PTA_).frobInner(this->densityA_->conjugate()).real();
 
    if(!this->isClosedShell && this->Ref_ != TCS){
      this->energyOneE += 
        (*this->aointegrals_->oneE_).frobInner(this->densityB_->real());
      this->energyTwoE += 
        0.5*(*this->PTB_).frobInner(this->densityB_->conjugate()).real();
    }
    // Add in the electric field component if they are non-zero
    std::array<double,3> null{{0,0,0}};
    if(this->elecField_ != null){
      int NB = this->nTCS_*this->nBasis_;
      int NBSq = NB*NB;
      int iBuf = 0;
      for(auto iXYZ = 0; iXYZ < 3; iXYZ++){
        ConstRealMap 
          mu(&this->aointegrals_->elecDipole_->storage()[iBuf],NB,NB);

        this->energyOneE += 
          this->elecField_[iXYZ] * mu.frobInner(this->densityA_->real());
        if(!this->isClosedShell && this->Ref_ != TCS)
        this->energyOneE += 
          this->elecField_[iXYZ] * mu.frobInner(this->densityB_->real());
        iBuf += NBSq;
      }
    }
 
 
    this->totalEnergy= this->energyOneE + this->energyTwoE + this->energyNuclei;
  }
#ifdef CQ_ENABLE_MPI
  MPI_Bcast(&this->totalEnergy,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
  MPI_Bcast(&this->energyOneE,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
  MPI_Bcast(&this->energyTwoE,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
#endif
//this->printEnergy();
};
*/

template<>
void SingleSlater<dcomplex>::getAlgebraicField(){ 
  this->algebraicField_      = "Complex";
  this->algebraicFieldShort_ = "\u2102";
}

template<>
void SingleSlater<dcomplex>::writeSCFFiles(){

  this->fileio_->alphaSCFDen->write(this->densityA_->data(),
    *(this->fileio_->complexType)
  );
  this->fileio_->alphaMO->write(this->moA_->data(),
    *(this->fileio_->complexType)
  );
  if(!this->isClosedShell && this->Ref_ != TCS){
    this->fileio_->betaSCFDen->write(this->densityB_->data(),
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
} // Namespace ChronusQ
