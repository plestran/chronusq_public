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
//---------------------
// initialize MOIntegrals
//---------------------
void MOIntegrals::iniMOIntegrals(   Molecule * molecule, BasisSet * basisSet, FileIO * fileio, 
                                    Controls * controls, AOIntegrals * aointegrals, SingleSlater<double> * singleSlater){
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


void MOIntegrals::getLocMO(){
  if(this->haveLocMO) return;

  if(this->Ref_ == SingleSlater<double>::TCS){
    this->locMOOcc_ = 
      std::unique_ptr<RealTensor2d>(new RealTensor2d(this->nTCS_*this->nBasis_,this->nO_));
    this->locMOVir_ = 
      std::unique_ptr<RealTensor2d>(new RealTensor2d(this->nTCS_*this->nBasis_,this->nV_));
  } else {
    this->locMOAOcc_ = 
      std::unique_ptr<RealTensor2d>(new RealTensor2d(this->nBasis_,this->nOA_));
    this->locMOAVir_ = 
      std::unique_ptr<RealTensor2d>(new RealTensor2d(this->nBasis_,this->nVA_));
    if(!this->singleSlater_->isClosedShell){
      this->locMOBOcc_ = 
        std::unique_ptr<RealTensor2d>(new RealTensor2d(this->nBasis_,this->nOB_));
      this->locMOBVir_ = 
        std::unique_ptr<RealTensor2d>(new RealTensor2d(this->nBasis_,this->nVB_));
    }
  }

  if(this->Ref_ == SingleSlater<double>::TCS){
    for(auto i = 0; i < this->nTCS_*this->nBasis_; i++)
    for(auto j = 0; j < this->nO_                ; j++){
      (*this->locMOOcc_)(i,j) = (*this->singleSlater_->moA())(i,j);
    }
    for(auto i = 0        ; i < this->nTCS_*this->nBasis_; i++)
    for(auto j = this->nO_; j < this->nTCS_*this->nBasis_; j++){
      (*this->locMOOcc_)(i,j - this->nO_) = (*this->singleSlater_->moA())(i,j);
    }
  } else {
    for(auto i = 0; i < this->nBasis_; i++)
    for(auto j = 0; j < this->nOA_   ; j++){
      (*this->locMOAOcc_)(i,j) = (*this->singleSlater_->moA())(i,j);
    }
    for(auto i = 0         ; i < this->nBasis_; i++)
    for(auto j = this->nOA_; j < this->nBasis_; j++){
      (*this->locMOAOcc_)(i,j - this->nOA_) = (*this->singleSlater_->moA())(i,j);
    }
    if(!this->singleSlater_->isClosedShell){
      for(auto i = 0; i < this->nBasis_; i++)
      for(auto j = 0; j < this->nOB_   ; j++){
        (*this->locMOBOcc_)(i,j) = (*this->singleSlater_->moB())(i,j);
      }
      for(auto i = 0         ; i < this->nBasis_; i++)
      for(auto j = this->nOB_; j < this->nBasis_; j++){
        (*this->locMOBOcc_)(i,j - this->nOB_) = (*this->singleSlater_->moB())(i,j);
      }
    }
  }

  this->haveLocMO = true;
}

void MOIntegrals::formIAJB(bool doDBar){
  if(this->haveMOiajb && (this->iajbIsDBar == doDBar)) return;
  else if(this->iajbIsDBar != doDBar) this->iajb_.reset();
  this->getLocMO();

  RealTensor4d IinlsA,    IinlsB,    IinlsTCS;
  RealTensor4d IialsAA,   IialsBB,   IialsTCS;
  RealTensor4d IiajsAAA,  IiajsAAB,  IiajsBBB,  IiajsTCS;
  RealTensor4d SiajbAAAA, SiajbAABB, SiajbBBBB, SiajbTCS;

  if(this->Ref_ == SingleSlater<double>::TCS){
    auto NTCSxNBASIS = this->nTCS_*this->nBasis_;
    this->iajb_ = 
      std::unique_ptr<RealTensor4d>(new 
        RealTensor4d(this->nO_,this->nV_,this->nO_,this->nV_));

    IinlsTCS = RealTensor4d(this->nO_,NTCSxNBASIS,NTCSxNBASIS,NTCSxNBASIS);
    IialsTCS = RealTensor4d(this->nO_,this->nV_,NTCSxNBASIS,NTCSxNBASIS);
    IiajsTCS = RealTensor4d(this->nO_,this->nV_,this->nO_,NTCSxNBASIS);
    if(doDBar) SiajbTCS = RealTensor4d(this->nO_,this->nV_,this->nO_,this->nV_);
  } else {
    this->iajbAAAA_ = 
      std::unique_ptr<RealTensor4d>(new 
        RealTensor4d(this->nOA_,this->nVA_,this->nOA_,this->nVA_));
    this->iajbAABB_ = 
      std::unique_ptr<RealTensor4d>(new 
        RealTensor4d(this->nOA_,this->nVA_,this->nOB_,this->nVB_));
    if(!this->singleSlater_->isClosedShell) {
      this->iajbBBBB_ = 
        std::unique_ptr<RealTensor4d>(new 
          RealTensor4d(this->nOB_,this->nVB_,this->nOB_,this->nVB_));
    }

    IinlsA   = RealTensor4d(this->nOA_,this->nBasis_,this->nBasis_,this->nBasis_);
    IialsAA  = RealTensor4d(this->nOA_,this->nVA_,this->nBasis_,this->nBasis_);
    IiajsAAA = RealTensor4d(this->nOA_,this->nVA_,this->nOA_,this->nBasis_);
    if(doDBar) {
      SiajbAAAA = RealTensor4d(this->nOA_,this->nVA_,this->nOA_,this->nVA_);
    }
    if(!this->singleSlater_->isClosedShell) {
      IinlsB   = RealTensor4d(this->nOB_,this->nBasis_,this->nBasis_,this->nBasis_);
      IialsBB  = RealTensor4d(this->nOB_,this->nVB_,this->nBasis_,this->nBasis_);
      IiajsAAB = RealTensor4d(this->nOA_,this->nVA_,this->nOB_,this->nBasis_);
      IiajsBBB = RealTensor4d(this->nOB_,this->nVB_,this->nOB_,this->nBasis_);
      if(doDBar) {
        SiajbAABB = RealTensor4d(this->nOA_,this->nVA_,this->nOB_,this->nVB_);
        SiajbBBBB = RealTensor4d(this->nOB_,this->nVB_,this->nOB_,this->nVB_);
      }
    }
  }

  enum{i,j,k,l,a,b,c,d,m,n,lm,s}; 

  if(this->Ref_ == SingleSlater<double>::TCS){
    contract(1.0,(*this->locMOOcc_),{m,i},(*this->aointegrals_->aoERI_),{m,n,lm,s},
             0.0,IinlsTCS,{i,n,lm,s});
    contract(1.0,(*this->locMOVir_),{n,a},IinlsTCS,{i,n,lm,s},0.0,IialsTCS,{i,a,lm,s});
    contract(1.0,(*this->locMOOcc_),{lm,j},IialsTCS,{i,a,lm,s},0.0,IiajsTCS,{i,a,j,s});
    if(doDBar)
      contract(1.0,(*this->locMOVir_),{s,b},IiajsTCS,{i,a,j,s},0.0,SiajbTCS,{i,a,j,b});
    else
      contract(1.0,(*this->locMOVir_),{s,b},IiajsTCS,{i,a,j,s},0.0,(*this->iajb_),{i,a,j,b});
  } else {
    contract(1.0,(*this->locMOAOcc_),{m,i},(*this->aointegrals_->aoERI_),{m,n,lm,s},
             0.0,IinlsA,{i,n,lm,s});
    contract(1.0,(*this->locMOAVir_),{n,a},IinlsA,{i,n,lm,s},0.0,IialsAA,{i,a,lm,s});
    contract(1.0,(*this->locMOAOcc_),{lm,j},IialsAA,{i,a,lm,s},0.0,IiajsAAA,{i,a,j,s});

    if(doDBar)
      contract(1.0,(*this->locMOAVir_),{s,b},IiajsAAA,{i,a,j,s},0.0,SiajbAAAA,{i,a,j,b});
    else
      contract(1.0,(*this->locMOAVir_),{s,b},IiajsAAA,{i,a,j,s},
               0.0,(*this->iajbAAAA_),{i,a,j,b});

    if(!this->singleSlater_->isClosedShell){
      contract(1.0,(*this->locMOBOcc_),{m,i},(*this->aointegrals_->aoERI_),{m,n,lm,s},
               0.0,IinlsB,{i,n,lm,s});
      contract(1.0,(*this->locMOBVir_),{n,a},IinlsB,{i,n,lm,s},0.0,IialsBB,{i,a,lm,s});
      contract(1.0,(*this->locMOBOcc_),{lm,j},IialsAA,{i,a,lm,s},0.0,IiajsAAB,{i,a,j,s});
      contract(1.0,(*this->locMOBOcc_),{lm,j},IialsBB,{i,a,lm,s},0.0,IiajsBBB,{i,a,j,s});

      contract(1.0,(*this->locMOBVir_),{s,b},IiajsAAB,{i,a,j,s},
               0.0,(*this->iajbAABB_),{i,a,j,b});
      if(doDBar){
        contract(1.0,(*this->locMOBVir_),{s,b},IiajsBBB,{i,a,j,s},0.0,SiajbBBBB,{i,a,j,b});
      } else {
        contract(1.0,(*this->locMOBVir_),{s,b},IiajsBBB,{i,a,j,s},
                 0.0,(*this->iajbBBBB_),{i,a,j,b});
      }
    } else {
      contract(1.0,(*this->locMOAVir_),{s,b},IiajsAAA,{i,a,j,s},
               0.0,(*this->iajbAABB_),{i,a,j,b});
    }
  }

  if(doDBar){
    if(this->Ref_ == SingleSlater<double>::TCS){
      for(auto i = 0; i < this->nO_; i++)
      for(auto a = 0; a < this->nV_; a++)
      for(auto j = 0; j < this->nO_; j++)
      for(auto b = 0; b < this->nV_; b++) 
        (*this->iajb_)(i,a,j,b) = SiajbTCS(i,a,j,b) - SiajbTCS(i,b,j,a);
    } else {
      for(auto i = 0; i < this->nOA_; i++)
      for(auto a = 0; a < this->nVA_; a++)
      for(auto j = 0; j < this->nOA_; j++)
      for(auto b = 0; b < this->nVA_; b++) 
        (*this->iajbAAAA_)(i,a,j,b) = SiajbAAAA(i,a,j,b) - SiajbAAAA(i,b,j,a);

      if(!this->singleSlater_->isClosedShell)
        for(auto i = 0; i < this->nOB_; i++)
        for(auto a = 0; a < this->nVB_; a++)
        for(auto j = 0; j < this->nOB_; j++)
        for(auto b = 0; b < this->nVB_; b++) 
          (*this->iajbBBBB_)(i,a,j,b) = SiajbBBBB(i,a,j,b) - SiajbBBBB(i,b,j,a);
    }
  }

  
  this->haveMOiajb = true;
  this->iajbIsDBar = doDBar;
}
