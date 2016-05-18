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
//APS
//  Function to read and match the Gaussian Guess from previous calculation
//  Note: It works as long you use the cartisian basis in both softwares
template<>
void SingleSlater<double>::matchord(){
  unsigned int ibase=0;
  unsigned int irow;
  for(auto i=0;i<this->nShell_;i++){
    int L=basisset_->shells(i).contr[0].l;
//  s functions order matchs
    if (L==0){ibase += 1;}
//  p functions order matchs
    else if (L==1){ibase+=3;} 
//  d functions in cartesian --> swap
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
//  f functions in cartesian --> swap
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
//  g functions in cartesian --> swap
    else if (L==4){
      irow = ibase-1;
      for(auto ishift=0;ishift<7;ishift++) {
        (*this->moA_).row((irow+(15-ishift))).swap(
            (*this->moA_).row(irow+(1+ishift))
          );
        if(this->Ref_ != RHF) 
          (*this->moB_).row((irow+(15-ishift))).swap(
              (*this->moB_).row(irow+(1+ishift))
          );
      }
      ibase+=15;
    }
//  h function in cartesian --> swap
    else if (L==5){
      irow = ibase-1;
      for(auto ishift=0;ishift<10;ishift++) {
        (*this->moA_).row((irow+(21-ishift))).swap(
            (*this->moA_).row(irow+(1+ishift))
          );
        if(this->Ref_ != RHF) 
          (*this->moB_).row((irow+(21-ishift))).swap(
              (*this->moB_).row(irow+(1+ishift))
            );
      }
      ibase+=21;
    }
//  I function in cartesian --> swap
    else if (L==6){
      irow = ibase-1;
      for(auto ishift=0;ishift<14;ishift++) {
        (*this->moA_).row((irow+(28-ishift))).swap(
            (*this->moA_).row(irow+(1+ishift))
          );
        if(this->Ref_ != RHF) 
          (*this->moB_).row((irow+(28-ishift))).swap(
              (*this->moB_).row(irow+(1+ishift))
            );
      }
      ibase+=28;
    }
//    K function in cartesian --> swap
    else if (L==7){
      irow = ibase-1;
      for(auto ishift=0;ishift<18;ishift++) {
        (*this->moA_).row((irow+(36-ishift))).swap(
            (*this->moA_).row(irow+(1+ishift))
          );
        if(this->Ref_ != RHF) 
          (*this->moB_).row((irow+(36-ishift))).swap(
              (*this->moB_).row(irow+(1+ishift))
            );
      }
      ibase+=36;
    }
    else {
      this->fileio_->out<<"L>7? Really? not yet implemented"<<endl;
    };
  };
  prettyPrint(this->fileio_->out,(*this->moA_),"APS PRINT GUESS SWAP");
};

template<>
void SingleSlater<double>::getAlgebraicField(){ 
  this->algebraicField_      = "Real";
  this->algebraicFieldShort_ = "\u211D";
}

template<>
void SingleSlater<double>::writeSCFFiles(){
  this->fileio_->alphaSCFDen->write(this->onePDMA_->data(),
      H5::PredType::NATIVE_DOUBLE);
  this->fileio_->alphaMO->write(this->moA_->data(),
      H5::PredType::NATIVE_DOUBLE);
  if(!this->isClosedShell && this->Ref_ != TCS){
    this->fileio_->betaSCFDen->write(this->onePDMB_->data(),
        H5::PredType::NATIVE_DOUBLE);
    this->fileio_->betaMO->write(this->moB_->data(),
        H5::PredType::NATIVE_DOUBLE);
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

template<>
void SingleSlater<double>::backTransformMOs(){
  if(this->nTCS_ == 1) {
    this->NBSqScratch_->noalias() = 
      (*this->aointegrals_->ortho1_) * (*this->moA_);
    (*this->moA_) = (*this->NBSqScratch_);

    if(!this->isClosedShell){
      this->NBSqScratch_->noalias() = 
        (*this->aointegrals_->ortho1_) * (*this->moB_);
      (*this->moB_) = (*this->NBSqScratch_);
    }
  } else {
    Eigen::Map<RealMatrix,0,Eigen::Stride<Dynamic,Dynamic> >
      MOA(this->moA_->data(),this->nBasis_,this->nTCS_*this->nBasis_,
          Eigen::Stride<Dynamic,Dynamic>(this->nTCS_*this->nBasis_,2));
    Eigen::Map<RealMatrix,0,Eigen::Stride<Dynamic,Dynamic> >
      MOB(this->moA_->data()+1,this->nBasis_,this->nTCS_*this->nBasis_,
          Eigen::Stride<Dynamic,Dynamic>(this->nTCS_*this->nBasis_,2));

    RealMap SCRATCH1(this->memManager_->malloc<double>(this->nBasis_*
          this->nBasis_*this->nTCS_),this->nBasis_,this->nTCS_*this->nBasis_);
    RealMap SCRATCH2(this->memManager_->malloc<double>(this->nBasis_*
          this->nBasis_*this->nTCS_),this->nBasis_,this->nTCS_*this->nBasis_);

    SCRATCH1 = MOA;
    SCRATCH2 = (*this->aointegrals_->ortho1_) * SCRATCH1;
    MOA = SCRATCH2; 

    SCRATCH1 = MOB;
    SCRATCH2 = (*this->aointegrals_->ortho1_) * SCRATCH1;
    MOB = SCRATCH2; 

    this->memManager_->free(SCRATCH1.data(),
        this->nBasis_*this->nBasis_*this->nTCS_);
    this->memManager_->free(SCRATCH2.data(),
        this->nBasis_*this->nBasis_*this->nTCS_);
    
  }
};

} // Namespace ChronusQ
