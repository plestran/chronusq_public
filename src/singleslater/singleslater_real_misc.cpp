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
SingleSlater<double>::SingleSlater(SingleSlater<double> * other){
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
    this->densityA_           = std::unique_ptr<RealMatrix>(
      new RealMatrix(*other->densityA_)
    );
    
    if(getRank() == 0) {
      this->fockA_              = std::unique_ptr<RealMatrix>(
        new RealMatrix(*other->fockA_)
      );
      this->moA_                = std::unique_ptr<RealMatrix>(
        new RealMatrix(*other->moA_)
      );
      this->PTA_                = std::unique_ptr<RealMatrix>(
        new RealMatrix(*other->PTA_)
      );
    }
    if(this->Ref_ != isClosedShell && this->Ref_ != TCS ){
      this->densityB_           = std::unique_ptr<RealMatrix>(
        new RealMatrix(*other->densityB_)
      );

      if(getRank() == 0) {
        this->fockB_              = std::unique_ptr<RealMatrix>(
          new RealMatrix(*other->fockB_)
        );
        this->moB_                = std::unique_ptr<RealMatrix>(
          new RealMatrix(*other->moB_)
        );
        this->PTB_                = std::unique_ptr<RealMatrix>(
          new RealMatrix(*other->PTB_)
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

//APS
//  Function to read and match the Gaussian Guess from previous calculation
//  Note: It works as long you use the cartisian basis in both softwares
  template<>
  void SingleSlater<double>::matchord(){
    unsigned int ibase=0;
    unsigned int irow;
    for(auto i=0;i<this->nShell_;i++){
          int L=basisset_->shells(i).contr[0].l;
//        cout << L << endl;
//         s functions order matchs
           if (L==0){ibase += 1;}
//         p functions order matchs
           else if (L==1){ibase+=3;} 
//           d functions in cartesian --> swap
           else if (L==2){
             irow = ibase-1;
             (*this->moA_).row((irow+2)).swap((*this->moA_).row(irow+4));
             (*this->moA_).row((irow+3)).swap((*this->moA_).row(irow+5));
             (*this->moA_).row((irow+5)).swap((*this->moA_).row(irow+6));
             if(this->Ref_ != RHF){
               (*this->moB_).row((irow+2)).swap((*this->moB_).row(irow+4));
               (*this->moB_).row((irow+3)).swap((*this->moB_).row(irow+5));
               (*this->moB_).row((irow+5)).swap((*this->moB_).row(irow+6));
             }
             ibase+=6;
             }
//           f functions in cartesian --> swap
           else if (L==3){
             irow = ibase-1;
             (*this->moA_).row((irow+2)).swap((*this->moA_).row(irow+5));
             (*this->moA_).row((irow+3)).swap((*this->moA_).row(irow+6));
             (*this->moA_).row((irow+5)).swap((*this->moA_).row(irow+10));
             (*this->moA_).row((irow+6)).swap((*this->moA_).row(irow+7));
             (*this->moA_).row((irow+7)).swap((*this->moA_).row(irow+10));
             (*this->moA_).row((irow+8)).swap((*this->moA_).row(irow+9));
             if(this->Ref_ != RHF) {
               (*this->moB_).row((irow+2)).swap((*this->moB_).row(irow+5));
               (*this->moB_).row((irow+3)).swap((*this->moB_).row(irow+6));
               (*this->moB_).row((irow+5)).swap((*this->moB_).row(irow+10));
               (*this->moB_).row((irow+6)).swap((*this->moB_).row(irow+7));
               (*this->moB_).row((irow+7)).swap((*this->moB_).row(irow+10));
               (*this->moB_).row((irow+8)).swap((*this->moB_).row(irow+9));
             }
             ibase+=10;
             }
//           g functions in cartesian --> swap
           else if (L==4){
             irow = ibase-1;
             for(auto ishift=0;ishift<7;ishift++) {
               (*this->moA_).row((irow+(15-ishift))).swap((*this->moA_).row(irow+(1+ishift)));
               if(this->Ref_ != RHF) (*this->moB_).row((irow+(15-ishift))).swap((*this->moB_).row(irow+(1+ishift)));
               }
             ibase+=15;
             }
//           h function in cartesian --> swap
           else if (L==5){
             irow = ibase-1;
             for(auto ishift=0;ishift<10;ishift++) {
               (*this->moA_).row((irow+(21-ishift))).swap((*this->moA_).row(irow+(1+ishift)));
               if(this->Ref_ != RHF) (*this->moB_).row((irow+(21-ishift))).swap((*this->moB_).row(irow+(1+ishift)));
               }
             ibase+=21;
             }
//           I function in cartesian --> swap
           else if (L==6){
             irow = ibase-1;
             for(auto ishift=0;ishift<14;ishift++) {
               (*this->moA_).row((irow+(28-ishift))).swap((*this->moA_).row(irow+(1+ishift)));
               if(this->Ref_ != RHF) (*this->moB_).row((irow+(28-ishift))).swap((*this->moB_).row(irow+(1+ishift)));
               }
             ibase+=28;
             }
//           K function in cartesian --> swap
           else if (L==7){
             irow = ibase-1;
             for(auto ishift=0;ishift<18;ishift++) {
               (*this->moA_).row((irow+(36-ishift))).swap((*this->moA_).row(irow+(1+ishift)));
               if(this->Ref_ != RHF) (*this->moB_).row((irow+(36-ishift))).swap((*this->moB_).row(irow+(1+ishift)));
               }
             ibase+=36;
             }
           else {this->fileio_->out<<"L>7? Really? not yet implemented"<<endl;};
       };
     prettyPrint(this->fileio_->out,(*this->moA_),"APS PRINT GUESS SWAP");
   };
//
/************************
 * Compute Total Energy *
 ************************/
template<>
void SingleSlater<double>::computeEnergy(){
  if(getRank() == 0) {
    double energyXC;
    this->energyOneE = 
      (*this->aointegrals_->oneE_).frobInner(this->densityA_->conjugate());
    this->energyTwoE = 
      0.5*(*this->PTA_).frobInner(this->densityA_->conjugate());
 
    if(!this->isClosedShell && this->Ref_ != TCS){
      this->energyOneE += 
        (*this->aointegrals_->oneE_).frobInner(this->densityB_->conjugate());
      this->energyTwoE += 
        0.5*(*this->PTB_).frobInner(this->densityB_->conjugate());
    }
    
    if(this->isDFT) this->energyTwoE += this->totalEx + this->totalEcorr;
      
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
          this->elecField_[iXYZ] * mu.frobInner(this->densityA_->conjugate());

        if(!this->isClosedShell && this->Ref_ != TCS)
          this->energyOneE += 
            this->elecField_[iXYZ] * mu.frobInner(this->densityB_->conjugate());
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

template<>
void SingleSlater<double>::getAlgebraicField(){ 
  this->algebraicField_      = "Real";
  this->algebraicFieldShort_ = "\u211D";
}

template<>
void SingleSlater<double>::writeSCFFiles(){
  this->fileio_->alphaSCFDen->write(this->densityA_->data(),H5::PredType::NATIVE_DOUBLE);
  this->fileio_->alphaMO->write(this->moA_->data(),H5::PredType::NATIVE_DOUBLE);
  if(!this->isClosedShell && this->Ref_ != TCS){
    this->fileio_->betaSCFDen->write(this->densityB_->data(),H5::PredType::NATIVE_DOUBLE);
    this->fileio_->betaMO->write(this->moB_->data(),H5::PredType::NATIVE_DOUBLE);
  }
}

template<>
void SingleSlater<double>::fixPhase(){
  if(!this->fixPhase_) return;
   RealMatrix::Index maxIndex;
   for(auto iCol = 0; iCol < this->moA_->cols(); iCol++){
     this->moA_->col(iCol).cwiseAbs().maxCoeff(&maxIndex);
     if(this->moA_->col(iCol)(maxIndex) < 0)
       this->moA_->col(iCol) *= -1;
   }
};

} // Namespace ChronusQ
