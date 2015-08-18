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
#include <mointegrals.h>
using ChronusQ::AOIntegrals;
using ChronusQ::Molecule;
using ChronusQ::BasisSet;
using ChronusQ::Controls;
using ChronusQ::FileIO;
using ChronusQ::MOIntegrals;
namespace ChronusQ{
//---------------------
// initialize MOIntegrals
//---------------------
template<>
void MOIntegrals<dcomplex>::iniMOIntegrals(   Molecule * molecule, BasisSet * basisSet, FileIO * fileio, 
                                    Controls * controls, AOIntegrals * aointegrals, SingleSlater<dcomplex> * singleSlater){
  this->molecule_ = molecule;
  this->basisSet_ = basisSet;
  this->fileio_   = fileio;
  this->controls_ = controls;
  this->aointegrals_ = aointegrals;
  this->singleSlater_ = singleSlater;

  this->haveMOiajb = false;
  this->haveMOijab = false;
  this->haveMOijka = false;
  this->haveMOijkl = false;
  this->haveMOiabc = false;
  this->haveMOabcd = false;
  this->iajbIsDBar = false;
  this->ijabIsDBar = false;
  this->ijkaIsDBar = false;
  this->ijklIsDBar = false;
  this->iabcIsDBar = false;
  this->abcdIsDBar = false;
  this->haveLocMO  = false;

  this->nBasis_ = singleSlater->nBasis();
  this->Ref_    = singleSlater->Ref();
  this->nTCS_   = singleSlater->nTCS();
  this->nOA_    = singleSlater->nOccA();
  this->nOB_    = singleSlater->nOccB();
  this->nVA_    = singleSlater->nVirA();
  this->nVB_    = singleSlater->nVirB();
  this->nO_     = this->nOA_ + this->nOB_;  
  this->nV_     = this->nVA_ + this->nVB_;
};

template<>
void MOIntegrals<dcomplex>::getLocMO(){
  if(this->haveLocMO) return;

  if(this->Ref_ == SingleSlater<dcomplex>::TCS){
    this->reLocMOOcc_ = 
      std::unique_ptr<RealTensor2d>(
        new RealTensor2d(this->nTCS_*this->nBasis_,this->nO_));
    this->reLocMOVir_ = 
      std::unique_ptr<RealTensor2d>(
        new RealTensor2d(this->nTCS_*this->nBasis_,this->nV_));
    this->imLocMOOcc_ = 
      std::unique_ptr<RealTensor2d>(
        new RealTensor2d(this->nTCS_*this->nBasis_,this->nO_));
    this->imLocMOVir_ = 
      std::unique_ptr<RealTensor2d>(
        new RealTensor2d(this->nTCS_*this->nBasis_,this->nV_));
  } else {
    this->reLocMOAOcc_ = 
      std::unique_ptr<RealTensor2d>(
        new RealTensor2d(this->nBasis_,this->nOA_));
    this->reLocMOAVir_ = 
      std::unique_ptr<RealTensor2d>(
        new RealTensor2d(this->nBasis_,this->nVA_));
    this->imLocMOAOcc_ = 
      std::unique_ptr<RealTensor2d>(
        new RealTensor2d(this->nBasis_,this->nOA_));
    this->imLocMOAVir_ = 
      std::unique_ptr<RealTensor2d>(
        new RealTensor2d(this->nBasis_,this->nVA_));
    if(!this->singleSlater_->isClosedShell){
      this->reLocMOBOcc_ = 
        std::unique_ptr<RealTensor2d>(
          new RealTensor2d(this->nBasis_,this->nOB_));
      this->reLocMOBVir_ = 
        std::unique_ptr<RealTensor2d>(
          new RealTensor2d(this->nBasis_,this->nVB_));
      this->imLocMOBOcc_ = 
        std::unique_ptr<RealTensor2d>(
          new RealTensor2d(this->nBasis_,this->nOB_));
      this->imLocMOBVir_ = 
        std::unique_ptr<RealTensor2d>(
          new RealTensor2d(this->nBasis_,this->nVB_));
    }
  }

  if(this->Ref_ == SingleSlater<double>::TCS){
    for(auto i = 0; i < this->nTCS_*this->nBasis_; i++)
    for(auto j = 0; j < this->nO_                ; j++){
      (*this->locMOOcc_)(i,j) = (*this->singleSlater_->moA())(i,j);
    }
    for(auto i = 0        ; i < this->nTCS_*this->nBasis_; i++)
    for(auto j = this->nO_; j < this->nTCS_*this->nBasis_; j++){
      (*this->locMOVir_)(i,j - this->nO_) = (*this->singleSlater_->moA())(i,j);
    }
  } else {
    for(auto i = 0; i < this->nBasis_; i++)
    for(auto j = 0; j < this->nOA_   ; j++){
      (*this->locMOAOcc_)(i,j) = (*this->singleSlater_->moA())(i,j);
    }
    for(auto i = 0         ; i < this->nBasis_; i++)
    for(auto j = this->nOA_; j < this->nBasis_; j++){
      (*this->locMOAVir_)(i,j - this->nOA_) = (*this->singleSlater_->moA())(i,j);
    }
    if(!this->singleSlater_->isClosedShell){
      for(auto i = 0; i < this->nBasis_; i++)
      for(auto j = 0; j < this->nOB_   ; j++){
        (*this->locMOBOcc_)(i,j) = (*this->singleSlater_->moB())(i,j);
      }
      for(auto i = 0         ; i < this->nBasis_; i++)
      for(auto j = this->nOB_; j < this->nBasis_; j++){
        (*this->locMOBVir_)(i,j - this->nOB_) = (*this->singleSlater_->moB())(i,j);
      }
    }
  }

  if(this->Ref_ == SingleSlater<dcomplex>::TCS){
    for(auto i = 0; i < this->nTCS_*this->nBasis_; i++)
    for(auto j = 0; j < this->nO_                ; j++){
      (*this->reLocMOOcc_)(i,j) = (*this->locMOOcc_)(i,j).real();
      (*this->imLocMOOcc_)(i,j) = (*this->locMOOcc_)(i,j).imag();
    }
    for(auto i = 0; i < this->nTCS_*this->nBasis_; i++)
    for(auto j = 0; j < this->nV_                ; j++){
      (*this->reLocMOVir_)(i,j) = (*this->locMOVir_)(i,j+this->nO_).real();
      (*this->imLocMOVir_)(i,j) = (*this->locMOVir_)(i,j+this->nO_).imag();
    }
  } else {
    for(auto i = 0; i < this->nTCS_*this->nBasis_ ; i++)
    for(auto j = 0; j < this->nOA_                ; j++){
      (*this->reLocMOAOcc_)(i,j) = (*this->locMOAOcc_)(i,j).real();
      (*this->imLocMOAOcc_)(i,j) = (*this->locMOAOcc_)(i,j).imag();
    }
    for(auto i = 0; i < this->nTCS_*this->nBasis_ ; i++)
    for(auto j = 0; j < this->nVA_                ; j++){
      (*this->reLocMOAVir_)(i,j) = (*this->locMOAVir_)(i,j+this->nOA_).real();
      (*this->imLocMOAVir_)(i,j) = (*this->locMOAVir_)(i,j+this->nOA_).imag();
    }
    if(!this->singleSlater_->isClosedShell){
      for(auto i = 0; i < this->nTCS_*this->nBasis_ ; i++)
      for(auto j = 0; j < this->nOB_                ; j++){
        (*this->reLocMOBOcc_)(i,j) = (*this->locMOBOcc_)(i,j).real();
        (*this->imLocMOBOcc_)(i,j) = (*this->locMOBOcc_)(i,j).imag();
      }
      for(auto i = 0; i < this->nTCS_*this->nBasis_ ; i++)
      for(auto j = 0; j < this->nVB_                ; j++){
        (*this->reLocMOBVir_)(i,j) = (*this->locMOBVir_)(i,j+this->nOB_).real();
        (*this->imLocMOBVir_)(i,j) = (*this->locMOBVir_)(i,j+this->nOB_).imag();
      }
    }
  }

  this->haveLocMO = true;
}

}//namespace ChronusQ
