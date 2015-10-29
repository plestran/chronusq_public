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
template<>
void MOIntegrals<double>::getLocMO(){
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

  this->haveLocMO = true;
}
/*
void MOIntegrals::formIAJB(bool doDBar){
  if(this->haveMOiajb && (this->iajbIsDBar == doDBar)) return;
  else if(this->iajbIsDBar != doDBar) this->iajb_.reset();
  cout << "HERE 1" << endl;
  this->getLocMO();
  cout << "HERE 1" << endl;

  RealTensor4d IinlsA,    IinlsB,    IinlsTCS;
  RealTensor4d IialsAA,   IialsBB,   IialsTCS;
  RealTensor4d IiajsAAA,  IiajsAAB,  IiajsBBB,  IiajsTCS;
  RealTensor4d SiajbAAAA, SiajbAABB, SiajbBBBB, SiajbTCS;
  cout << "HERE 1" << endl;

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
  cout << "HERE 1" << endl;
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

  cout << "HERE 1" << endl;
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
  cout << "HERE 2" << endl;
    contract(1.0,(*this->locMOAOcc_),{m,i},(*this->aointegrals_->aoERI_),{m,n,lm,s},
             0.0,IinlsA,{i,n,lm,s});
  cout << "HERE 2" << endl;
    contract(1.0,(*this->locMOAVir_),{n,a},IinlsA,{i,n,lm,s},0.0,IialsAA,{i,a,lm,s});
  cout << "HERE 2" << endl;
    contract(1.0,(*this->locMOAOcc_),{lm,j},IialsAA,{i,a,lm,s},0.0,IiajsAAA,{i,a,j,s});
  cout << "HERE 2" << endl;

    if(doDBar)
      contract(1.0,(*this->locMOAVir_),{s,b},IiajsAAA,{i,a,j,s},0.0,SiajbAAAA,{i,a,j,b});
    else
      contract(1.0,(*this->locMOAVir_),{s,b},IiajsAAA,{i,a,j,s},
               0.0,(*this->iajbAAAA_),{i,a,j,b});
  cout << "HERE 2" << endl;

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
  cout << "HERE 2" << endl;
      contract(1.0,(*this->locMOAVir_),{s,b},IiajsAAA,{i,a,j,s},
               0.0,(*this->iajbAABB_),{i,a,j,b});
  cout << "HERE 2" << endl;
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
*/
template<>
void MOIntegrals<double>::testLocMO(){
  this->getLocMO();
  RealTensor2d P,S,C,I;

  P = RealTensor2d(this->nTCS_*this->nBasis_,this->nTCS_*this->nBasis_);
  C = RealTensor2d(this->nTCS_*this->nBasis_,this->nTCS_*this->nBasis_);
  S = RealTensor2d(this->nTCS_*this->nBasis_,this->nTCS_*this->nBasis_);
  I = RealTensor2d(this->nTCS_*this->nBasis_,this->nTCS_*this->nBasis_);

  for(auto i = 0; i < this->nTCS_*this->nBasis_; i++)
  for(auto j = 0; j < this->nTCS_*this->nBasis_; j++){
    S(i,j) = (*this->aointegrals_->overlap_)(i,j);
  }
  
  for(auto i = 0; i < this->nTCS_*this->nBasis_; i++)
  for(auto j = 0; j < this->nOA_               ; j++){
    C(i,j) = (*this->locMOAOcc_)(i,j);
  }
  for(auto i = 0         ; i < this->nTCS_*this->nBasis_; i++)
  for(auto j = this->nOA_; j < this->nTCS_*this->nBasis_; j++){
    C(i,j) = (*this->locMOAVir_)(i,j-this->nOA_);
  }

  enum{mm,nn,ii,jj};
  contract(2.0,(*this->locMOAOcc_),{mm,ii},(*this->locMOAOcc_),{nn,ii},0.0,P,{mm,nn}); 

  cout << "Density Diff" << endl;
  for(auto i = 0; i < this->nTCS_*this->nBasis_; i++){
  for(auto j = 0; j < this->nTCS_*this->nBasis_; j++){
  cout << P(i,j) - (*this->singleSlater_->densityA())(i,j) << "   "; 
  }
  cout << endl;
  }

  contract(1.0,C,{mm,ii},S,{mm,nn},0.0,P,{ii,nn}); 
  contract(1.0,P,{ii,nn},C,{nn,jj},0.0,I,{ii,jj}); 
  cout << endl << endl << "C**H S C" << endl;
  for(auto i = 0; i < this->nTCS_*this->nBasis_; i++){
  for(auto j = 0; j < this->nTCS_*this->nBasis_; j++){
  cout << I(i,j)  << "   "; 
  }
  cout << endl;
  }
}

template<>
void MOIntegrals<double>::formIAJB(bool doDBar){
  if(this->haveMOiajb && (this->iajbIsDBar == doDBar)) return;
  else if(this->iajbIsDBar != doDBar) {
    if(this->Ref_ == SingleSlater<double>::TCS)
      this->iajb_.reset();
    else {
      this->iajbAAAA_.reset();
      this->iajbAABB_.reset();
      if(!this->singleSlater_->isClosedShell) 
        this->iajbBBBB_.reset();
    }
  }

  this->getLocMO();


  if(this->Ref_ == SingleSlater<double>::TCS){
    this->iajb_ = std::unique_ptr<RealTensor4d>(
      new RealTensor4d(this->nO_,this->nV_,this->nO_,this->nV_));
  } else {
    this->iajbAAAA_ = std::unique_ptr<RealTensor4d>(
      new RealTensor4d(this->nOA_,this->nVA_,this->nOA_,this->nVA_));
    this->iajbAABB_ = std::unique_ptr<RealTensor4d>(
      new RealTensor4d(this->nOA_,this->nVA_,this->nOB_,this->nVB_));
    if(!this->singleSlater_->isClosedShell)
      this->iajbBBBB_ = std::unique_ptr<RealTensor4d>(
        new RealTensor4d(this->nOB_,this->nVB_,this->nOB_,this->nVB_));
  }

  RealTensor4d Iinls , Iials  , Iiajs   ;
  RealTensor4d Siajb ;
  RealTensor4d IinlsA, IialsAA, IiajsAAA;  
  RealTensor4d IinlsB, IialsBB, IiajsBBB;
  RealTensor4d IiajsAAB, SiajbAAAA, SiajbBBBB;
  auto NTCSxNBASIS = this->nTCS_*this->nBasis_;

  if(this->Ref_ == SingleSlater<double>::TCS){
    Iinls = RealTensor4d(this->nO_,NTCSxNBASIS,NTCSxNBASIS,NTCSxNBASIS);
    Iials = RealTensor4d(this->nO_,this->nV_,NTCSxNBASIS,NTCSxNBASIS);
    Iiajs = RealTensor4d(this->nO_,this->nV_,this->nO_,NTCSxNBASIS);
    if(doDBar)
      Siajb  = RealTensor4d(this->nO_,this->nV_,this->nO_,this->nV_);
  } else {
    IinlsA    = RealTensor4d(this->nOA_,this->nBasis_,this->nBasis_,this->nBasis_);
    IialsAA   = RealTensor4d(this->nOA_,this->nVA_,this->nBasis_,this->nBasis_);
    IiajsAAA  = RealTensor4d(this->nOA_,this->nVA_,this->nOA_,this->nBasis_);
    if(!this->singleSlater_->isClosedShell){
      IinlsB     = RealTensor4d(this->nOB_,this->nBasis_,this->nBasis_,this->nBasis_);
      IialsBB    = RealTensor4d(this->nOB_,this->nVB_,this->nBasis_,this->nBasis_);
      IiajsBBB   = RealTensor4d(this->nOB_,this->nVB_,this->nOB_,this->nBasis_);
      IiajsAAB   = RealTensor4d(this->nOA_,this->nVA_,this->nOB_,this->nBasis_);
      if(doDBar){
        SiajbAAAA  = RealTensor4d(this->nOA_,this->nVA_,this->nOA_,this->nVA_);
        SiajbBBBB  = RealTensor4d(this->nOB_,this->nVB_,this->nOB_,this->nVB_);
      }
    }
  }

  enum{mu,nu,lm,sg,i,a,j,b};

  if(this->Ref_ == SingleSlater<double>::TCS){

    // First Quarter Transformation  (i ν | λ κ)
    contract(1.0,(*this->locMOOcc_),{mu,i},(*this->aointegrals_->aoERI_),{mu,nu,lm,sg},
             0.0,Iinls,{i,nu,lm,sg});
 
    // First Half Transformation     (i a | λ κ)
    contract(1.0,(*this->locMOVir_),{nu,a},Iinls,{i,nu,lm,sg},0.0,Iials,{i,a,lm,sg});
 
    // Third Quarter Transformation  (i a | j κ)
    contract(1.0,(*this->locMOOcc_),{lm,j},Iials,{i,a,lm,sg},0.0,Iiajs,{i,a,j,sg});

  } else {

    // First Quarter Transformation Alpha             (i ν | λ κ) [A]
    contract(1.0,(*this->locMOAOcc_),{mu,i},(*this->aointegrals_->aoERI_),{mu,nu,lm,sg},
             0.0,IinlsA,{i,nu,lm,sg});
 
    // First Half Transformation Alpha-Alpha          (i a | λ κ) [AA]
    contract(1.0,(*this->locMOAVir_),{nu,a},IinlsA,{i,nu,lm,sg},0.0,IialsAA,{i,a,lm,sg});
 
    // Third Quarter Transformation Alpha-Alpha-Alpha (i a | j κ) [AA|A]
    contract(1.0,(*this->locMOAOcc_),{lm,j},IialsAA,{i,a,lm,sg},0.0,IiajsAAA,{i,a,j,sg});
 
    /******************************/
    /* ONLY BUILD IF CLOSED SHELL */
    /******************************/
    if(!this->singleSlater_->isClosedShell){
 
      // First Quarter Transformation Beta            (i ν | λ κ) [B]
      contract(1.0,(*this->locMOBOcc_),{mu,i},(*this->aointegrals_->aoERI_),{mu,nu,lm,sg},
               0.0,IinlsB,{i,nu,lm,sg});
 
      // First Half Transformation Beta-Beta          (i a | λ κ) [BB]
      contract(1.0,(*this->locMOBVir_),{nu,a},IinlsB,{i,nu,lm,sg},0.0,IialsBB,{i,a,lm,sg});
 
      // Third Quarter Transformation Beta-Beta-Beta    (i a | j κ) [BB|B]
      contract(1.0,(*this->locMOBOcc_),{lm,j},IialsBB,{i,a,lm,sg},0.0,IiajsBBB,{i,a,j,sg});
 
      // Third Quarter Transformation Alpha-Alpha-Beta  (i a | j κ) [AA|B]
      contract(1.0,(*this->locMOBOcc_),{lm,j},IialsAA,{i,a,lm,sg},0.0,IiajsAAB,{i,a,j,sg});
 
    }
  }

  /******************************/
  /*  IF ONLY DOING SINGLE BAR  */
  /******************************/
  if(!doDBar){
    /*
     * If only doing single bar integrals, populate the pure spin
     * storage with the single bar MO integrals
     */ 
 
    if(this->Ref_ == SingleSlater<double>::TCS){
      // Last Quarter Transformation (i a | j b) 
      contract(1.0,(*this->locMOVir_),{sg,b},Iiajs,{i,a,j,sg},
               0.0,(*this->iajb_),{i,a,j,b});
    } else {

      // Last Quarter Transformation Alpha-Alpha-Alpha-Alpha (i a | j b) [AA|AA]
      contract(1.0,(*this->locMOAVir_),{sg,b},IiajsAAA,{i,a,j,sg},
               0.0,(*this->iajbAAAA_),{i,a,j,b});
     
      /******************************/
      /* ONLY BUILD IF CLOSED SHELL */
      /******************************/
      // Last Quarter Transformation Beta-Beta-Beta-Beta (i a | j b) [BB|BB]
      if(!this->singleSlater_->isClosedShell)
        contract(1.0,(*this->locMOBVir_),{sg,b},IiajsBBB,{i,a,j,sg},
                 0.0,(*this->iajbBBBB_),{i,a,j,b});
    }

  /****************************************/
  /*  IF DOING DOUBLE AND OPEN-SHELL BAR  */
  /****************************************/
  } else if(!this->singleSlater_->isClosedShell || this->Ref_ == SingleSlater<double>::TCS) {
    /*
     * If doing double bar integrals with an open-shell or two-component reference, 
     * populate previously allocated intermetiates with the single bar MO integrals
     */ 

    if(this->Ref_ == SingleSlater<double>::TCS)
      // Last Quarter Transformation (i a | j b)
      contract(1.0,(*this->locMOVir_),{sg,b},Iiajs,{i,a,j,sg},0.0,Siajb,{i,a,j,b});
    else {
      // Last Quarter Transformation Alpha-Alpha-Alpha-Alpha (i a | j b) [AA|AA]
      contract(1.0,(*this->locMOAVir_),{sg,b},IiajsAAA,{i,a,j,sg},0.0,SiajbAAAA,{i,a,j,b});
      // Last Quarter Transformation Beta-Beta-Beta-Beta (i a | j b) [BB|BB]
      contract(1.0,(*this->locMOBVir_),{sg,b},IiajsBBB,{i,a,j,sg},0.0,SiajbBBBB,{i,a,j,b});
    }
  }

  if(this->Ref_ != SingleSlater<double>::TCS){
    /*
     * Regardless of whether or not we are doing double-bar integrals, populate
     * the mixed spin storage with single-bar integrals (for RHF and UHF, the spin
     * integration removes the exchange term)
     */ 
    if(!this->singleSlater_->isClosedShell){
      // (UHF) Last Quarter Transformation Alpha-Alpha-Beta-Beta (i a | j b) [AA|BB]
      contract(1.0,(*this->locMOBVir_),{sg,b},IiajsAAB,{i,a,j,sg},
               0.0,(*this->iajbAABB_),{i,a,j,b});
    } else {
      // (RHF) Last Quarter Transformation Alpha-Alpha-Beta-Beta (i a | j b) [AA|BB]
      contract(1.0,(*this->locMOAVir_),{sg,b},IiajsAAA,{i,a,j,sg},
               0.0,(*this->iajbAABB_),{i,a,j,b});
    }
  }

  /*
   * Build double-bar integrals from single-bar integrals is requested
   *
   * (i a || j b) = (i a | j b) - (i b | j a)
   *
   */
  if(doDBar){
    if(this->Ref_ == SingleSlater<double>::TCS)
      for(auto i = 0; i < this->nO_; i++)
      for(auto a = 0; a < this->nV_; a++)
      for(auto j = 0; j < this->nO_; j++)
      for(auto b = 0; b < this->nV_; b++)
        (*this->iajb_)(i,a,j,b) = Siajb(i,a,j,b) - Siajb(i,b,j,a);
    else {
      for(auto i = 0; i < this->nOA_; i++)
      for(auto a = 0; a < this->nVA_; a++)
      for(auto j = 0; j < this->nOA_; j++)
      for(auto b = 0; b < this->nVA_; b++){
        /*
         * In the case of RHF, we can reuse the "mixed-spin" storage to build the
         * double bar integrals, as the mixed-spin storage is simply the single
         * bar integrals. This does not work for UHF
         */ 
        if(this->singleSlater_->isClosedShell)
          (*this->iajbAAAA_)(i,a,j,b) =
            (*this->iajbAABB_)(i,a,j,b)-(*this->iajbAABB_)(i,b,j,a);
        else {
          (*this->iajbAAAA_)(i,a,j,b) = SiajbAAAA(i,a,j,b) - SiajbAAAA(i,b,j,a);
          (*this->iajbBBBB_)(i,a,j,b) = SiajbBBBB(i,a,j,b) - SiajbBBBB(i,b,j,a);
        }
      }
      if(!this->singleSlater_->isClosedShell)
        for(auto i = 0; i < this->nOB_; i++)
        for(auto a = 0; a < this->nVB_; a++)
        for(auto j = 0; j < this->nOB_; j++)
        for(auto b = 0; b < this->nVB_; b++)
          (*this->iajbBBBB_)(i,a,j,b) = SiajbBBBB(i,a,j,b) - SiajbBBBB(i,b,j,a);
    }
  }
/*
    for(auto i = 0; i < this->nOA_; i++)
    for(auto a = 0; a < this->nVA_; a++)
    for(auto j = 0; j < this->nOB_; j++)
    for(auto b = 0; b < this->nVB_; b++){
  cout << "(" << i << " " << a + this->nOA_ << " | " << j << " " << this->nOA_+b << ") " <<(*this->iajbAABB_)(i,a,j,b) << endl;
    }
*/
}

template<>
void MOIntegrals<double>::formABCD(bool doDBar){
  if(this->haveMOabcd && (this->abcdIsDBar == doDBar)) return;
  else if(this->abcdIsDBar != doDBar) {
    if(this->Ref_ == SingleSlater<double>::TCS)
      this->abcd_.reset();
    else {
      this->abcdAAAA_.reset();
      this->abcdAABB_.reset();
      if(!this->singleSlater_->isClosedShell) 
        this->abcdBBBB_.reset();
    }
  }

  this->getLocMO();

  if(this->Ref_ == SingleSlater<double>::TCS){
    this->abcd_ = std::unique_ptr<RealTensor4d>(
      new RealTensor4d(this->nV_,this->nV_,this->nV_,this->nV_));
  } else {
    this->abcdAAAA_ = std::unique_ptr<RealTensor4d>(
      new RealTensor4d(this->nVA_,this->nVA_,this->nVA_,this->nVA_));
    this->abcdAABB_ = std::unique_ptr<RealTensor4d>(
      new RealTensor4d(this->nVA_,this->nVA_,this->nVB_,this->nVB_));
    if(!this->singleSlater_->isClosedShell)
      this->abcdBBBB_ = std::unique_ptr<RealTensor4d>(
        new RealTensor4d(this->nVB_,this->nVB_,this->nVB_,this->nVB_));
  }

  RealTensor4d Ianls , Iabls  , Iabcs   ;
  RealTensor4d Sabcd ;
  RealTensor4d IanlsA, IablsAA, IabcsAAA;  
  RealTensor4d IanlsB, IablsBB, IabcsBBB;
  RealTensor4d IabcsAAB, SabcdAAAA, SabcdBBBB;
  auto NTCSxNBASIS = this->nTCS_*this->nBasis_;

  if(this->Ref_ == SingleSlater<double>::TCS){
    Ianls = RealTensor4d(this->nV_,NTCSxNBASIS,NTCSxNBASIS,NTCSxNBASIS);
    Iabls = RealTensor4d(this->nV_,this->nV_,NTCSxNBASIS,NTCSxNBASIS);
    Iabcs = RealTensor4d(this->nV_,this->nV_,this->nV_,NTCSxNBASIS);
    if(doDBar)
      Sabcd  = RealTensor4d(this->nV_,this->nV_,this->nV_,this->nV_);
  } else {
    IanlsA    = RealTensor4d(this->nVA_,this->nBasis_,this->nBasis_,this->nBasis_);
    IablsAA   = RealTensor4d(this->nVA_,this->nVA_,this->nBasis_,this->nBasis_);
    IabcsAAA  = RealTensor4d(this->nVA_,this->nVA_,this->nVA_,this->nBasis_);
    if(!this->singleSlater_->isClosedShell){
      IanlsB     = RealTensor4d(this->nVB_,this->nBasis_,this->nBasis_,this->nBasis_);
      IablsBB    = RealTensor4d(this->nVB_,this->nVB_,this->nBasis_,this->nBasis_);
      IabcsBBB   = RealTensor4d(this->nVB_,this->nVB_,this->nVB_,this->nBasis_);
      IabcsAAB   = RealTensor4d(this->nVA_,this->nVA_,this->nVB_,this->nBasis_);
      if(doDBar){
        SabcdAAAA  = RealTensor4d(this->nVA_,this->nVA_,this->nVA_,this->nVA_);
        SabcdBBBB  = RealTensor4d(this->nVB_,this->nVB_,this->nVB_,this->nVB_);
      }
    }
  }

  enum{mu,nu,lm,sg,a,b,c,d};

  if(this->Ref_ == SingleSlater<double>::TCS){

    // First Quarter Transformation  (a ν | λ κ)
    contract(1.0,(*this->locMOVir_),{mu,a},(*this->aointegrals_->aoERI_),{mu,nu,lm,sg},
             0.0,Ianls,{a,nu,lm,sg});
 
    // First Half Transformation     (a b | λ κ)
    contract(1.0,(*this->locMOVir_),{nu,b},Ianls,{a,nu,lm,sg},0.0,Iabls,{a,b,lm,sg});
 
    // Third Quarter Transformation  (a b | c κ)
    contract(1.0,(*this->locMOVir_),{lm,c},Iabls,{a,b,lm,sg},0.0,Iabcs,{a,b,c,sg});

  } else {

    // First Quarter Transformation Alpha             (a ν | λ κ) [A]
    contract(1.0,(*this->locMOAVir_),{mu,a},(*this->aointegrals_->aoERI_),{mu,nu,lm,sg},
             0.0,IanlsA,{a,nu,lm,sg});
 
    // First Half Transformation Alpha-Alpha          (a b | λ κ) [AA]
    contract(1.0,(*this->locMOAVir_),{nu,b},IanlsA,{a,nu,lm,sg},0.0,IablsAA,{a,b,lm,sg});
 
    // Third Quarter Transformation Alpha-Alpha-Alpha (a b | c κ) [AA|A]
    contract(1.0,(*this->locMOAVir_),{lm,c},IablsAA,{a,b,lm,sg},0.0,IabcsAAA,{a,b,c,sg});
 
    /******************************/
    /* ONLY BUILD IF CLOSED SHELL */
    /******************************/
    if(!this->singleSlater_->isClosedShell){
 
      // First Quarter Transformation Beta            (a ν | λ κ) [B]
      contract(1.0,(*this->locMOBVir_),{mu,a},(*this->aointegrals_->aoERI_),{mu,nu,lm,sg},
               0.0,IanlsB,{a,nu,lm,sg});
 
      // First Half Transformation Beta-Beta          (a b | λ κ) [BB]
      contract(1.0,(*this->locMOBVir_),{nu,b},IanlsB,{a,nu,lm,sg},0.0,IablsBB,{a,b,lm,sg});
 
      // Third Quarter Transformation Beta-Beta-Beta    (a b | c κ) [BB|B]
      contract(1.0,(*this->locMOBVir_),{lm,c},IablsBB,{a,b,lm,sg},0.0,IabcsBBB,{a,b,c,sg});
 
      // Third Quarter Transformation Alpha-Alpha-Beta  (a b | c κ) [AA|B]
      contract(1.0,(*this->locMOBVir_),{lm,c},IablsAA,{a,b,lm,sg},0.0,IabcsAAB,{a,b,c,sg});
 
    }
  }

  /******************************/
  /*  IF ONLY DOING SINGLE BAR  */
  /******************************/
  if(!doDBar){
    /*
     * If only doing single bar integrals, populate the pure spin
     * storage with the single bar MO integrals
     */ 
 
    if(this->Ref_ == SingleSlater<double>::TCS){
      // Last Quarter Transformation (a b | c d) 
      contract(1.0,(*this->locMOVir_),{sg,d},Iabcs,{a,b,c,sg},
               0.0,(*this->abcd_),{a,b,c,d});
    } else {

      // Last Quarter Transformation Alpha-Alpha-Alpha-Alpha (a b | c d) [AA|AA]
      contract(1.0,(*this->locMOAVir_),{sg,d},IabcsAAA,{a,b,c,sg},
               0.0,(*this->abcdAAAA_),{a,b,c,d});
     
      /******************************/
      /* ONLY BUILD IF CLOSED SHELL */
      /******************************/
      // Last Quarter Transformation Beta-Beta-Beta-Beta (a b | c d) [BB|BB]
      if(!this->singleSlater_->isClosedShell)
        contract(1.0,(*this->locMOBVir_),{sg,d},IabcsBBB,{a,b,c,sg},
                 0.0,(*this->abcdBBBB_),{a,b,c,d});
    }

  /****************************************/
  /*  IF DOING DOUBLE AND OPEN-SHELL BAR  */
  /****************************************/
  } else if(!this->singleSlater_->isClosedShell || this->Ref_ == SingleSlater<double>::TCS) {
    /*
     * If doing double bar integrals with an open-shell or two-component reference, 
     * populate previously allocated intermetiates with the single bar MO integrals
     */ 

    if(this->Ref_ == SingleSlater<double>::TCS)
      // Last Quarter Transformation (a b | c d)
      contract(1.0,(*this->locMOVir_),{sg,d},Iabcs,{a,b,c,sg},0.0,Sabcd,{a,b,c,d});
    else {
      // Last Quarter Transformation Alpha-Alpha-Alpha-Alpha (a b | c d) [AA|AA]
      contract(1.0,(*this->locMOAVir_),{sg,d},IabcsAAA,{a,b,c,sg},0.0,SabcdAAAA,{a,b,c,d});
      // Last Quarter Transformation Beta-Beta-Beta-Beta (a b | c d) [BB|BB]
      contract(1.0,(*this->locMOBVir_),{sg,d},IabcsBBB,{a,b,c,sg},0.0,SabcdBBBB,{a,b,c,d});
    }
  }

  if(this->Ref_ != SingleSlater<double>::TCS){
    /*
     * Regardless of whether or not we are doing double-bar integrals, populate
     * the mixed spin storage with single-bar integrals (for RHF and UHF, the spin
     * integration removes the exchange term)
     */ 
    if(!this->singleSlater_->isClosedShell){
      // (UHF) Last Quarter Transformation Alpha-Alpha-Beta-Beta (a b | c d) [AA|BB]
      contract(1.0,(*this->locMOBVir_),{sg,d},IabcsAAB,{a,b,c,sg},
               0.0,(*this->abcdAABB_),{a,b,c,d});
    } else {
      // (RHF) Last Quarter Transformation Alpha-Alpha-Beta-Beta (a b | c d) [AA|BB]
      contract(1.0,(*this->locMOAVir_),{sg,d},IabcsAAA,{a,b,c,sg},
               0.0,(*this->abcdAABB_),{a,b,c,d});
    }
  }

  /*
   * Build double-bar integrals from single-bar integrals is requested
   *
   * (a b || c d) = (a b | c d) - (a d | c b)
   *
   */
  if(doDBar){
    if(this->Ref_ == SingleSlater<double>::TCS)
      for(auto a = 0; a < this->nV_; a++)
      for(auto b = 0; b < this->nV_; b++)
      for(auto c = 0; c < this->nV_; c++)
      for(auto d = 0; d < this->nV_; d++)
        (*this->abcd_)(a,b,c,d) = Sabcd(a,b,c,d) - Sabcd(a,d,c,b);
    else {
      for(auto a = 0; a < this->nVA_; a++)
      for(auto b = 0; b < this->nVA_; b++)
      for(auto c = 0; c < this->nVA_; c++)
      for(auto d = 0; d < this->nVA_; d++){
        /*
         * In the case of RHF, we can reuse the "mixed-spin" storage to build the
         * double bar integrals, as the mixed-spin storage is simply the single
         * bar integrals. This does not work for UHF
         */ 
        if(this->singleSlater_->isClosedShell)
          (*this->abcdAAAA_)(a,b,c,d) =
            (*this->abcdAABB_)(a,b,c,d)-(*this->abcdAABB_)(a,d,c,b);
        else {
          (*this->abcdAAAA_)(a,b,c,d) = SabcdAAAA(a,b,c,d) - SabcdAAAA(a,d,c,b);
          (*this->abcdBBBB_)(a,b,c,d) = SabcdBBBB(a,b,c,d) - SabcdBBBB(a,d,c,b);
        }
      }
      if(!this->singleSlater_->isClosedShell)
      for(auto a = 0; a < this->nVB_; a++)
      for(auto b = 0; b < this->nVB_; b++)
      for(auto c = 0; c < this->nVB_; c++)
      for(auto d = 0; d < this->nVB_; d++)
          (*this->abcdBBBB_)(a,b,c,d) = SabcdBBBB(a,b,c,d) - SabcdBBBB(a,d,c,b);
    }
  }
    
/*
    for(auto i = 0; i < this->nOA_; i++)
    for(auto a = 0; a < this->nVA_; a++)
    for(auto j = 0; j < this->nOB_; j++)
    for(auto b = 0; b < this->nVB_; b++){
  cout << "(" << i << " " << a + this->nOA_ << " | " << j << " " << this->nOA_+b << ") " <<(*this->iajbAABB_)(i,a,j,b) << endl;
    }
*/
}

template<>
void MOIntegrals<double>::formIJKL(bool doDBar){
  if(this->haveMOijkl && (this->ijklIsDBar == doDBar)) return;
  else if(this->ijklIsDBar != doDBar) {
    if(this->Ref_ == SingleSlater<double>::TCS)
      this->ijkl_.reset();
    else {
      this->ijklAAAA_.reset();
      this->ijklAABB_.reset();
      if(!this->singleSlater_->isClosedShell) 
        this->ijklBBBB_.reset();
    }
  }

  this->getLocMO();

  if(this->Ref_ == SingleSlater<double>::TCS){
    this->ijkl_ = std::unique_ptr<RealTensor4d>(
      new RealTensor4d(this->nO_,this->nO_,this->nO_,this->nO_));
  } else {
    this->ijklAAAA_ = std::unique_ptr<RealTensor4d>(
      new RealTensor4d(this->nOA_,this->nOA_,this->nOA_,this->nOA_));
    this->ijklAABB_ = std::unique_ptr<RealTensor4d>(
      new RealTensor4d(this->nOA_,this->nOA_,this->nOB_,this->nOB_));
    if(!this->singleSlater_->isClosedShell)
      this->ijklBBBB_ = std::unique_ptr<RealTensor4d>(
        new RealTensor4d(this->nOB_,this->nOB_,this->nOB_,this->nOB_));
  }

  RealTensor4d Iinls , Iijls  , Iijks   ;
  RealTensor4d Sijkl ;
  RealTensor4d IinlsA, IijlsAA, IijksAAA;  
  RealTensor4d IinlsB, IijlsBB, IijksBBB;
  RealTensor4d IijksAAB, SijklAAAA, SijklBBBB;
  auto NTCSxNBASIS = this->nTCS_*this->nBasis_;

  if(this->Ref_ == SingleSlater<double>::TCS){
    Iinls = RealTensor4d(this->nO_,NTCSxNBASIS,NTCSxNBASIS,NTCSxNBASIS);
    Iijls = RealTensor4d(this->nO_,this->nO_,NTCSxNBASIS,NTCSxNBASIS);
    Iijks = RealTensor4d(this->nO_,this->nO_,this->nO_,NTCSxNBASIS);
    if(doDBar)
      Sijkl  = RealTensor4d(this->nO_,this->nO_,this->nO_,this->nO_);
  } else {
    IinlsA    = RealTensor4d(this->nOA_,this->nBasis_,this->nBasis_,this->nBasis_);
    IijlsAA   = RealTensor4d(this->nOA_,this->nOA_,this->nBasis_,this->nBasis_);
    IijksAAA  = RealTensor4d(this->nOA_,this->nOA_,this->nOA_,this->nBasis_);
    if(!this->singleSlater_->isClosedShell){
      IinlsB     = RealTensor4d(this->nOB_,this->nBasis_,this->nBasis_,this->nBasis_);
      IijlsBB    = RealTensor4d(this->nOB_,this->nOB_,this->nBasis_,this->nBasis_);
      IijksBBB   = RealTensor4d(this->nOB_,this->nOB_,this->nOB_,this->nBasis_);
      IijksAAB   = RealTensor4d(this->nOA_,this->nOA_,this->nOB_,this->nBasis_);
      if(doDBar){
        SijklAAAA  = RealTensor4d(this->nOA_,this->nOA_,this->nOA_,this->nOA_);
        SijklBBBB  = RealTensor4d(this->nOB_,this->nOB_,this->nOB_,this->nOB_);
      }
    }
  }

  enum{mu,nu,lm,sg,i,j,k,l};

  if(this->Ref_ == SingleSlater<double>::TCS){

    // First Quarter Transformation  (a ν | λ κ)
    contract(1.0,(*this->locMOOcc_),{mu,i},(*this->aointegrals_->aoERI_),{mu,nu,lm,sg},
             0.0,Iinls,{i,nu,lm,sg});
 
    // First Half Transformation     (a b | λ κ)
    contract(1.0,(*this->locMOOcc_),{nu,j},Iinls,{i,nu,lm,sg},0.0,Iijls,{i,j,lm,sg});
 
    // Third Quarter Transformation  (a b | c κ)
    contract(1.0,(*this->locMOOcc_),{lm,k},Iijls,{i,j,lm,sg},0.0,Iijks,{i,j,k,sg});

  } else {

    // First Quarter Transformation Alpha             (a ν | λ κ) [A]
    contract(1.0,(*this->locMOAOcc_),{mu,i},(*this->aointegrals_->aoERI_),{mu,nu,lm,sg},
             0.0,IinlsA,{i,nu,lm,sg});
 
    // First Half Transformation Alpha-Alpha          (a b | λ κ) [AA]
    contract(1.0,(*this->locMOAOcc_),{nu,j},IinlsA,{i,nu,lm,sg},0.0,IijlsAA,{i,j,lm,sg});
 
    // Third Quarter Transformation Alpha-Alpha-Alpha (a b | c κ) [AA|A]
    contract(1.0,(*this->locMOAOcc_),{lm,k},IijlsAA,{i,j,lm,sg},0.0,IijksAAA,{i,j,k,sg});
 
    /******************************/
    /* ONLY BUILD IF CLOSED SHELL */
    /******************************/
    if(!this->singleSlater_->isClosedShell){
 
      // First Quarter Transformation Beta            (a ν | λ κ) [B]
      contract(1.0,(*this->locMOBOcc_),{mu,i},(*this->aointegrals_->aoERI_),{mu,nu,lm,sg},
               0.0,IinlsB,{i,nu,lm,sg});
 
      // First Half Transformation Beta-Beta          (a b | λ κ) [BB]
      contract(1.0,(*this->locMOBOcc_),{nu,j},IinlsB,{i,nu,lm,sg},0.0,IijlsBB,{i,j,lm,sg});
 
      // Third Quarter Transformation Beta-Beta-Beta    (a b | c κ) [BB|B]
      contract(1.0,(*this->locMOBOcc_),{lm,k},IijlsBB,{i,j,lm,sg},0.0,IijksBBB,{i,j,k,sg});
 
      // Third Quarter Transformation Alpha-Alpha-Beta  (a b | c κ) [AA|B]
      contract(1.0,(*this->locMOBOcc_),{lm,k},IijlsAA,{i,j,lm,sg},0.0,IijksAAB,{i,j,k,sg});
 
    }
  }

  /******************************/
  /*  IF ONLY DOING SINGLE BAR  */
  /******************************/
  if(!doDBar){
    /*
     * If only doing single bar integrals, populate the pure spin
     * storage with the single bar MO integrals
     */ 
 
    if(this->Ref_ == SingleSlater<double>::TCS){
      // Last Quarter Transformation (a b | c d) 
      contract(1.0,(*this->locMOOcc_),{sg,l},Iijks,{i,j,k,sg},
               0.0,(*this->ijkl_),{i,j,k,l});
    } else {

      // Last Quarter Transformation Alpha-Alpha-Alpha-Alpha (a b | c d) [AA|AA]
      contract(1.0,(*this->locMOAOcc_),{sg,l},IijksAAA,{i,j,k,sg},
               0.0,(*this->ijklAAAA_),{i,j,k,l});
     
      /******************************/
      /* ONLY BUILD IF CLOSED SHELL */
      /******************************/
      // Last Quarter Transformation Beta-Beta-Beta-Beta (a b | c d) [BB|BB]
      if(!this->singleSlater_->isClosedShell)
        contract(1.0,(*this->locMOBOcc_),{sg,l},IijksBBB,{i,j,k,sg},
                 0.0,(*this->ijklBBBB_),{i,j,k,l});
    }

  /****************************************/
  /*  IF DOING DOUBLE AND OPEN-SHELL BAR  */
  /****************************************/
  } else if(!this->singleSlater_->isClosedShell || this->Ref_ == SingleSlater<double>::TCS) {
    /*
     * If doing double bar integrals with an open-shell or two-component reference, 
     * populate previously allocated intermetiates with the single bar MO integrals
     */ 

    if(this->Ref_ == SingleSlater<double>::TCS)
      // Last Quarter Transformation (a b | c d)
      contract(1.0,(*this->locMOOcc_),{sg,l},Iijks,{i,j,k,sg},0.0,Sijkl,{i,j,k,l});
    else {
      // Last Quarter Transformation Alpha-Alpha-Alpha-Alpha (a b | c d) [AA|AA]
      contract(1.0,(*this->locMOAOcc_),{sg,l},IijksAAA,{i,j,k,sg},0.0,SijklAAAA,{i,j,k,l});
      // Last Quarter Transformation Beta-Beta-Beta-Beta (a b | c d) [BB|BB]
      contract(1.0,(*this->locMOBOcc_),{sg,l},IijksBBB,{i,j,k,sg},0.0,SijklBBBB,{i,j,k,l});
    }
  }

  if(this->Ref_ != SingleSlater<double>::TCS){
    /*
     * Regardless of whether or not we are doing double-bar integrals, populate
     * the mixed spin storage with single-bar integrals (for RHF and UHF, the spin
     * integration removes the exchange term)
     */ 
    if(!this->singleSlater_->isClosedShell){
      // (UHF) Last Quarter Transformation Alpha-Alpha-Beta-Beta (a b | c d) [AA|BB]
      contract(1.0,(*this->locMOBOcc_),{sg,l},IijksAAB,{i,j,k,sg},
               0.0,(*this->ijklAABB_),{i,j,k,l});
    } else {
      // (RHF) Last Quarter Transformation Alpha-Alpha-Beta-Beta (a b | c d) [AA|BB]
      contract(1.0,(*this->locMOAOcc_),{sg,l},IijksAAA,{i,j,k,sg},
               0.0,(*this->ijklAABB_),{i,j,k,l});
    }
  }

  /*
   * Build double-bar integrals from single-bar integrals is requested
   *
   * (a b || c d) = (a b | c d) - (a d | c b)
   *
   */
  if(doDBar){
    if(this->Ref_ == SingleSlater<double>::TCS)
      for(auto i = 0; i < this->nO_; i++)
      for(auto j = 0; j < this->nO_; j++)
      for(auto k = 0; k < this->nO_; k++)
      for(auto l = 0; l < this->nO_; l++)
        (*this->ijkl_)(i,j,k,l) = Sijkl(i,j,k,l) - Sijkl(i,l,k,j);
    else {
      for(auto i = 0; i < this->nOA_; i++)
      for(auto j = 0; j < this->nOA_; j++)
      for(auto k = 0; k < this->nOA_; k++)
      for(auto l = 0; l < this->nOA_; l++){
        /*
         * In the case of RHF, we can reuse the "mixed-spin" storage to build the
         * double bar integrals, as the mixed-spin storage is simply the single
         * bar integrals. This does not work for UHF
         */ 
        if(this->singleSlater_->isClosedShell)
          (*this->ijklAAAA_)(i,j,k,l) =
            (*this->ijklAABB_)(i,j,k,l)-(*this->ijklAABB_)(i,l,k,j);
        else {
          (*this->ijklAAAA_)(i,j,k,l) = SijklAAAA(i,j,k,l) - SijklAAAA(i,l,k,j);
          (*this->ijklBBBB_)(i,j,k,l) = SijklBBBB(i,j,k,l) - SijklBBBB(i,l,k,j);
        }
      }
      if(!this->singleSlater_->isClosedShell)
      for(auto i = 0; i < this->nOB_; i++)
      for(auto j = 0; j < this->nOB_; j++)
      for(auto k = 0; k < this->nOB_; k++)
      for(auto l = 0; l < this->nOB_; l++)
          (*this->ijklBBBB_)(i,j,k,l) = SijklBBBB(i,j,k,l) - SijklBBBB(i,l,k,j);
    }
  }
}

template<>
void MOIntegrals<double>::formIJAB(bool doDBar){
  if(this->haveMOijab && (this->ijabIsDBar == doDBar)) return;
  else if(this->ijabIsDBar != doDBar) {
    if(this->Ref_ == SingleSlater<double>::TCS)
      this->ijab_.reset();
    else {
      this->ijabAAAA_.reset();
      this->ijabAABB_.reset();
      if(!this->singleSlater_->isClosedShell) 
        this->ijabBBBB_.reset();
    }
  }

  this->getLocMO();


  if(this->Ref_ == SingleSlater<double>::TCS){
    this->ijab_ = std::unique_ptr<RealTensor4d>(
      new RealTensor4d(this->nO_,this->nO_,this->nV_,this->nV_));
  } else {
    this->ijabAAAA_ = std::unique_ptr<RealTensor4d>(
      new RealTensor4d(this->nOA_,this->nOA_,this->nVA_,this->nVA_));
    this->ijabAABB_ = std::unique_ptr<RealTensor4d>(
      new RealTensor4d(this->nOA_,this->nOA_,this->nVB_,this->nVB_));
    if(!this->singleSlater_->isClosedShell)
      this->ijabBBBB_ = std::unique_ptr<RealTensor4d>(
        new RealTensor4d(this->nOB_,this->nOB_,this->nVB_,this->nVB_));
  }

  RealTensor4d Iinls , Iijls  , Iijas   ;
  RealTensor4d Sijab ;
  RealTensor4d IinlsA, IijlsAA, IijasAAA;  
  RealTensor4d IinlsB, IijlsBB, IijasBBB;
  RealTensor4d IijasAAB, SijabAAAA, SijabBBBB;
  auto NTCSxNBASIS = this->nTCS_*this->nBasis_;

  if(this->Ref_ == SingleSlater<double>::TCS){
    Iinls = RealTensor4d(this->nO_,NTCSxNBASIS,NTCSxNBASIS,NTCSxNBASIS);
    Iijls = RealTensor4d(this->nO_,this->nO_,NTCSxNBASIS,NTCSxNBASIS);
    Iijas = RealTensor4d(this->nO_,this->nO_,this->nV_,NTCSxNBASIS);
    if(doDBar)
      Sijab  = RealTensor4d(this->nO_,this->nO_,this->nV_,this->nV_);
  } else {
    IinlsA    = RealTensor4d(this->nOA_,this->nBasis_,this->nBasis_,this->nBasis_);
    IijlsAA   = RealTensor4d(this->nOA_,this->nOA_,this->nBasis_,this->nBasis_);
    IijasAAA  = RealTensor4d(this->nOA_,this->nOA_,this->nVA_,this->nBasis_);
    if(!this->singleSlater_->isClosedShell){
      IinlsB     = RealTensor4d(this->nOB_,this->nBasis_,this->nBasis_,this->nBasis_);
      IijlsBB    = RealTensor4d(this->nOB_,this->nOB_,this->nBasis_,this->nBasis_);
      IijasBBB   = RealTensor4d(this->nOB_,this->nOB_,this->nVB_,this->nBasis_);
      IijasAAB   = RealTensor4d(this->nOA_,this->nOA_,this->nVB_,this->nBasis_);
      if(doDBar){
        SijabAAAA  = RealTensor4d(this->nOA_,this->nOA_,this->nVA_,this->nVA_);
        SijabBBBB  = RealTensor4d(this->nOB_,this->nOB_,this->nVB_,this->nVB_);
      }
    }
  }

  enum{mu,nu,lm,sg,i,a,j,b};

  if(this->Ref_ == SingleSlater<double>::TCS){

    // First Quarter Transformation  (i ν | λ κ)
    contract(1.0,(*this->locMOOcc_),{mu,i},(*this->aointegrals_->aoERI_),{mu,nu,lm,sg},
             0.0,Iinls,{i,nu,lm,sg});
 
    // First Half Transformation     (i j | λ κ)
    contract(1.0,(*this->locMOOcc_),{nu,j},Iinls,{i,nu,lm,sg},0.0,Iijls,{i,j,lm,sg});
 
    // Third Quarter Transformation  (i j | a κ)
    contract(1.0,(*this->locMOVir_),{lm,a},Iijls,{i,j,lm,sg},0.0,Iijas,{i,j,a,sg});

  } else {

    // First Quarter Transformation Alpha             (i ν | λ κ) [A]
    contract(1.0,(*this->locMOAOcc_),{mu,i},(*this->aointegrals_->aoERI_),{mu,nu,lm,sg},
             0.0,IinlsA,{i,nu,lm,sg});
 
    // First Half Transformation Alpha-Alpha          (i j | λ κ) [AA]
    contract(1.0,(*this->locMOAOcc_),{nu,j},IinlsA,{i,nu,lm,sg},0.0,IijlsAA,{i,j,lm,sg});
 
    // Third Quarter Transformation Alpha-Alpha-Alpha (i j | a κ) [AA|A]
    contract(1.0,(*this->locMOAVir_),{lm,a},IijlsAA,{i,j,lm,sg},0.0,IijasAAA,{i,j,a,sg});
 
    /******************************/
    /* ONLY BUILD IF CLOSED SHELL */
    /******************************/
    if(!this->singleSlater_->isClosedShell){
 
      // First Quarter Transformation Beta            (i ν | λ κ) [B]
      contract(1.0,(*this->locMOBOcc_),{mu,i},(*this->aointegrals_->aoERI_),{mu,nu,lm,sg},
               0.0,IinlsB,{i,nu,lm,sg});
 
      // First Half Transformation Beta-Beta          (i j | λ κ) [BB]
      contract(1.0,(*this->locMOBOcc_),{nu,j},IinlsB,{i,nu,lm,sg},0.0,IijlsBB,{i,j,lm,sg});
 
      // Third Quarter Transformation Beta-Beta-Beta    (i j | a κ) [BB|B]
      contract(1.0,(*this->locMOBVir_),{lm,a},IijlsBB,{i,j,lm,sg},0.0,IijasBBB,{i,j,a,sg});
 
      // Third Quarter Transformation Alpha-Alpha-Beta  (i j | a κ) [AA|B]
      contract(1.0,(*this->locMOBVir_),{lm,a},IijlsAA,{i,j,lm,sg},0.0,IijasAAB,{i,j,a,sg});
 
    }
  }

  /******************************/
  /*  IF ONLY DOING SINGLE BAR  */
  /******************************/
  if(!doDBar){
    /*
     * If only doing single bar integrals, populate the pure spin
     * storage with the single bar MO integrals
     */ 
 
    if(this->Ref_ == SingleSlater<double>::TCS){
      // Last Quarter Transformation (i j | a b) 
      contract(1.0,(*this->locMOVir_),{sg,b},Iijas,{i,j,a,sg},
               0.0,(*this->ijab_),{i,j,a,b});
    } else {

      // Last Quarter Transformation Alpha-Alpha-Alpha-Alpha (i j | a b) [AA|AA]
      contract(1.0,(*this->locMOAVir_),{sg,b},IijasAAA,{i,j,a,sg},
               0.0,(*this->ijabAAAA_),{i,j,a,b});
     
      /******************************/
      /* ONLY BUILD IF CLOSED SHELL */
      /******************************/
      // Last Quarter Transformation Beta-Beta-Beta-Beta (i j | a b) [BB|BB]
      if(!this->singleSlater_->isClosedShell)
        contract(1.0,(*this->locMOBVir_),{sg,b},IijasBBB,{i,j,a,sg},
                 0.0,(*this->ijabBBBB_),{i,j,a,b});
    }

  /****************************************/
  /*  IF DOING DOUBLE AND OPEN-SHELL BAR  */
  /****************************************/
  } else if(!this->singleSlater_->isClosedShell || this->Ref_ == SingleSlater<double>::TCS) {
    /*
     * If doing double bar integrals with an open-shell or two-component reference, 
     * populate previously allocated intermetiates with the single bar MO integrals
     */ 

    if(this->Ref_ == SingleSlater<double>::TCS)
      // Last Quarter Transformation (i j | a b)
      contract(1.0,(*this->locMOVir_),{sg,b},Iijas,{i,j,a,sg},0.0,Sijab,{i,j,a,b});
    else {
      // Last Quarter Transformation Alpha-Alpha-Alpha-Alpha (i j | a b) [AA|AA]
      contract(1.0,(*this->locMOAVir_),{sg,b},IijasAAA,{i,j,a,sg},0.0,SijabAAAA,{i,j,a,b});
      // Last Quarter Transformation Beta-Beta-Beta-Beta (i j | a b) [BB|BB]
      contract(1.0,(*this->locMOBVir_),{sg,b},IijasBBB,{i,j,a,sg},0.0,SijabBBBB,{i,j,a,b});
    }
  }

  if(this->Ref_ != SingleSlater<double>::TCS){
    /*
     * Regardless of whether or not we are doing double-bar integrals, populate
     * the mixed spin storage with single-bar integrals (for RHF and UHF, the spin
     * integration removes the exchange term)
     */ 
    if(!this->singleSlater_->isClosedShell){
      // (UHF) Last Quarter Transformation Alpha-Alpha-Beta-Beta (i j | a b) [AA|BB]
      contract(1.0,(*this->locMOBVir_),{sg,b},IijasAAB,{i,j,a,sg},
               0.0,(*this->ijabAABB_),{i,j,a,b});
    } else {
      // (RHF) Last Quarter Transformation Alpha-Alpha-Beta-Beta (i j | a b) [AA|BB]
      contract(1.0,(*this->locMOAVir_),{sg,b},IijasAAA,{i,j,a,sg},
               0.0,(*this->ijabAABB_),{i,j,a,b});
    }
  }
  if(doDBar){
    /// THIS IS WRONG --- DO NOT USE
    if(this->Ref_ == SingleSlater<double>::TCS)
      for(auto i = 0; i < this->nO_; i++)
      for(auto j = 0; j < this->nO_; j++)
      for(auto a = 0; a < this->nV_; a++)
      for(auto b = 0; b < this->nV_; b++)
        (*this->ijab_)(i,j,a,b) = Sijab(i,j,a,b) - Sijab(i,b,a,j);
    else {
      for(auto i = 0; i < this->nOA_; i++)
      for(auto j = 0; j < this->nOA_; j++)
      for(auto a = 0; a < this->nVA_; a++)
      for(auto b = 0; b < this->nVA_; b++){
        /*
         * In the case of RHF, we can reuse the "mixed-spin" storage to build the
         * double bar integrals, as the mixed-spin storage is simply the single
         * bar integrals. This does not work for UHF
         */ 
        if(this->singleSlater_->isClosedShell)
          (*this->ijabAAAA_)(i,j,a,b) =
            (*this->ijabAABB_)(i,j,a,b)-(*this->ijabAABB_)(i,b,a,j);
        else {
          (*this->ijabAAAA_)(i,j,a,b) = SijabAAAA(i,j,a,b) - SijabAAAA(i,b,a,j);
          (*this->ijabBBBB_)(i,j,a,b) = SijabBBBB(i,j,a,b) - SijabBBBB(i,b,a,j);
        }
      }
      if(!this->singleSlater_->isClosedShell)
        for(auto i = 0; i < this->nOB_; i++)
        for(auto j = 0; j < this->nOB_; j++)
        for(auto a = 0; a < this->nVB_; a++)
        for(auto b = 0; b < this->nVB_; b++)
          (*this->ijabBBBB_)(i,j,a,b) = SijabBBBB(i,j,a,b) - SijabBBBB(i,b,a,j);
    }
  }
}

template<>
//DO IABC
void MOIntegrals<double>::formIABC(bool doDBar){
    /// THIS IS WRONG --- DO NOT USE
  if(this->haveMOiabc && (this->iabcIsDBar == doDBar)) return;
  else if(this->iabcIsDBar != doDBar) {
    if(this->Ref_ == SingleSlater<double>::TCS)
      this->iabc_.reset();
    else {
      this->iabcAAAA_.reset();
      this->iabcAABB_.reset();
      if(!this->singleSlater_->isClosedShell) 
        this->iabcBBBB_.reset();
    }
  }

  this->getLocMO();


  if(this->Ref_ == SingleSlater<double>::TCS){
    this->iabc_ = std::unique_ptr<RealTensor4d>(
      new RealTensor4d(this->nO_,this->nV_,this->nV_,this->nV_));
  } else {
    this->iabcAAAA_ = std::unique_ptr<RealTensor4d>(
      new RealTensor4d(this->nOA_,this->nVA_,this->nVA_,this->nVA_));
    this->iabcAABB_ = std::unique_ptr<RealTensor4d>(
      new RealTensor4d(this->nOA_,this->nVA_,this->nVB_,this->nVB_));
    if(!this->singleSlater_->isClosedShell)
      this->iabcBBBB_ = std::unique_ptr<RealTensor4d>(
        new RealTensor4d(this->nOB_,this->nVB_,this->nVB_,this->nVB_));
  }


  RealTensor4d Iinls , Iials  , Iiabs   ;
  RealTensor4d Siabc ;
  RealTensor4d IinlsA, IialsAA, IiabsAAA;  
  RealTensor4d IinlsB, IialsBB, IiabsBBB;
  RealTensor4d IiabsAAB, SiabcAAAA, SiabcBBBB;
  auto NTCSxNBASIS = this->nTCS_*this->nBasis_;

  if(this->Ref_ == SingleSlater<double>::TCS){
    Iinls = RealTensor4d(this->nO_,NTCSxNBASIS,NTCSxNBASIS,NTCSxNBASIS);
    Iials = RealTensor4d(this->nO_,this->nV_,NTCSxNBASIS,NTCSxNBASIS);
    Iiabs = RealTensor4d(this->nO_,this->nV_,this->nV_,NTCSxNBASIS);
    if(doDBar)
      Siabc  = RealTensor4d(this->nO_,this->nV_,this->nV_,this->nV_);
  } else {
    IinlsA    = RealTensor4d(this->nOA_,this->nBasis_,this->nBasis_,this->nBasis_);
    IialsAA   = RealTensor4d(this->nOA_,this->nVA_,this->nBasis_,this->nBasis_);
    IiabsAAA  = RealTensor4d(this->nOA_,this->nVA_,this->nVA_,this->nBasis_);
    if(!this->singleSlater_->isClosedShell){
      IinlsB     = RealTensor4d(this->nOB_,this->nBasis_,this->nBasis_,this->nBasis_);
      IialsBB    = RealTensor4d(this->nOB_,this->nVB_,this->nBasis_,this->nBasis_);
      IiabsBBB   = RealTensor4d(this->nOB_,this->nVB_,this->nVB_,this->nBasis_);
      IiabsAAB   = RealTensor4d(this->nOA_,this->nVA_,this->nVB_,this->nBasis_);
      if(doDBar){
        SiabcAAAA  = RealTensor4d(this->nOA_,this->nVA_,this->nVA_,this->nVA_);
        SiabcBBBB  = RealTensor4d(this->nOB_,this->nVB_,this->nVB_,this->nVB_);
      }
    }
  }

  enum{mu,nu,lm,sg,i,b,a,c};

  if(this->Ref_ == SingleSlater<double>::TCS){

    // First Quarter Transformation  (i ν | λ κ)
    contract(1.0,(*this->locMOOcc_),{mu,i},(*this->aointegrals_->aoERI_),{mu,nu,lm,sg},
             0.0,Iinls,{i,nu,lm,sg});
 
    // First Half Transformation     (i a | λ κ)
    contract(1.0,(*this->locMOVir_),{nu,a},Iinls,{i,nu,lm,sg},0.0,Iials,{i,a,lm,sg});
 
    // Third Quarter Transformation  (i a | b κ)
    contract(1.0,(*this->locMOVir_),{lm,b},Iials,{i,a,lm,sg},0.0,Iiabs,{i,a,b,sg});

  } else {

    // First Quarter Transformation Alpha             (i ν | λ κ) [A]
    contract(1.0,(*this->locMOAOcc_),{mu,i},(*this->aointegrals_->aoERI_),{mu,nu,lm,sg},
             0.0,IinlsA,{i,nu,lm,sg});
 
    // First Half Transformation Alpha-Alpha          (i a | λ κ) [AA]
    contract(1.0,(*this->locMOAVir_),{nu,a},IinlsA,{i,nu,lm,sg},0.0,IialsAA,{i,a,lm,sg});
 
    // Third Quarter Transformation Alpha-Alpha-Alpha (i a | b κ) [AA|A]
    contract(1.0,(*this->locMOAVir_),{lm,b},IialsAA,{i,a,lm,sg},0.0,IiabsAAA,{i,a,b,sg});
 
    /******************************/
    /* ONLY BUILD IF CLOSED SHELL */
    /******************************/
    if(!this->singleSlater_->isClosedShell){
 
      // First Quarter Transformation Beta            (i ν | λ κ) [B]
      contract(1.0,(*this->locMOBOcc_),{mu,i},(*this->aointegrals_->aoERI_),{mu,nu,lm,sg},
               0.0,IinlsB,{i,nu,lm,sg});
 
      // First Half Transformation Beta-Beta          (i a | λ κ) [BB]
      contract(1.0,(*this->locMOBVir_),{nu,a},IinlsB,{i,nu,lm,sg},0.0,IialsBB,{i,a,lm,sg});
 
      // Third Quarter Transformation Beta-Beta-Beta    (i a | b κ) [BB|B]
      contract(1.0,(*this->locMOBVir_),{lm,b},IialsBB,{i,a,lm,sg},0.0,IiabsBBB,{i,a,b,sg});
 
      // Third Quarter Transformation Alpha-Alpha-Beta  (i j | a κ) [AA|B]
      contract(1.0,(*this->locMOBVir_),{lm,b},IialsAA,{i,a,lm,sg},0.0,IiabsAAB,{i,a,b,sg});
 
    }
  }

  /******************************/
  /*  IF ONLY DOING SINGLE BAR  */
  /******************************/
  if(!doDBar){
    /*
     * If only doing single bar integrals, populate the pure spin
     * storage with the single bar MO integrals
     */ 
 
    if(this->Ref_ == SingleSlater<double>::TCS){
      // Last Quarter Transformation (i a | b c) 
      contract(1.0,(*this->locMOVir_),{sg,c},Iiabs,{i,a,b,sg},
               0.0,(*this->iabc_),{i,a,b,c});
    } else {

      // Last Quarter Transformation Alpha-Alpha-Alpha-Alpha (i a | b c) [AA|AA]
      contract(1.0,(*this->locMOAVir_),{sg,c},IiabsAAA,{i,a,b,sg},
               0.0,(*this->iabcAAAA_),{i,a,b,c});
     
      /******************************/
      /* ONLY BUILD IF CLOSED SHELL */
      /******************************/
      // Last Quarter Transformation Beta-Beta-Beta-Beta (i a | b c) [BB|BB]
      if(!this->singleSlater_->isClosedShell)
        contract(1.0,(*this->locMOBVir_),{sg,c},IiabsBBB,{i,a,b,sg},
                 0.0,(*this->iabcBBBB_),{i,a,b,c});
    }

  /****************************************/
  /*  IF DOING DOUBLE AND OPEN-SHELL BAR  */
  /****************************************/
  } else if(!this->singleSlater_->isClosedShell || this->Ref_ == SingleSlater<double>::TCS) {
    /*
     * If doing double bar integrals with an open-shell or two-component reference, 
     * populate previously allocated intermetiates with the single bar MO integrals
     */ 

    if(this->Ref_ == SingleSlater<double>::TCS)
      // Last Quarter Transformation (i a | b c)
      contract(1.0,(*this->locMOVir_),{sg,c},Iiabs,{i,a,b,sg},0.0,Siabc,{i,a,b,c});
    else {
      // Last Quarter Transformation Alpha-Alpha-Alpha-Alpha (i a | b c) [AA|AA]
      contract(1.0,(*this->locMOAVir_),{sg,c},IiabsAAA,{i,a,b,sg},0.0,SiabcAAAA,{i,a,b,c});
      // Last Quarter Transformation Beta-Beta-Beta-Beta (i a | b c) [BB|BB]
      contract(1.0,(*this->locMOBVir_),{sg,c},IiabsBBB,{i,a,b,sg},0.0,SiabcBBBB,{i,a,b,c});
    }
  }

  if(this->Ref_ != SingleSlater<double>::TCS){
    /*
     * Regardless of whether or not we are doing double-bar integrals, populate
     * the mixed spin storage with single-bar integrals (for RHF and UHF, the spin
     * integration removes the exchange term)
     */ 
    if(!this->singleSlater_->isClosedShell){
      // (UHF) Last Quarter Transformation Alpha-Alpha-Beta-Beta (i a | b c) [AA|BB]
      contract(1.0,(*this->locMOBVir_),{sg,c},IiabsAAB,{i,a,b,sg},
               0.0,(*this->iabcAABB_),{i,a,b,c});
    } else {
      // (RHF) Last Quarter Transformation Alpha-Alpha-Beta-Beta (i a | b c) [AA|BB]
      contract(1.0,(*this->locMOAVir_),{sg,c},IiabsAAA,{i,a,b,sg},
               0.0,(*this->iabcAABB_),{i,a,b,c});
    }
  }
  if(doDBar){
    if(this->Ref_ == SingleSlater<double>::TCS)
      for(auto i = 0; i < this->nO_; i++)
      for(auto a = 0; a < this->nV_; a++)
      for(auto b = 0; b < this->nV_; b++)
      for(auto c = 0; c < this->nV_; c++)
        (*this->iabc_)(i,a,b,c) = Siabc(i,a,b,c) - Siabc(i,c,b,a);
    else {
      for(auto i = 0; i < this->nOA_; i++)
      for(auto a = 0; a < this->nVA_; a++)
      for(auto b = 0; b < this->nVA_; b++)
      for(auto c = 0; c < this->nVA_; c++){
        /*
         * In the case of RHF, we can reuse the "mixed-spin" storage to build the
         * double bar integrals, as the mixed-spin storage is simply the single
         * bar integrals. This does not work for UHF
         */ 
        if(this->singleSlater_->isClosedShell)
          (*this->iabcAAAA_)(i,a,b,c) =
            (*this->iabcAABB_)(i,a,b,c)-(*this->iabcAABB_)(i,c,b,a);
        else {
          (*this->iabcAAAA_)(i,a,b,c) = SiabcAAAA(i,a,b,c) - SiabcAAAA(i,c,b,a);
          (*this->iabcBBBB_)(i,a,b,c) = SiabcBBBB(i,a,b,c) - SiabcBBBB(i,c,b,a);
        }
      }
      if(!this->singleSlater_->isClosedShell)
        for(auto i = 0; i < this->nOB_; i++)
        for(auto a = 0; a < this->nVB_; a++)
        for(auto b = 0; b < this->nVB_; b++)
        for(auto c = 0; c < this->nVB_; c++)
          (*this->iabcBBBB_)(i,a,b,c) = SiabcBBBB(i,a,b,c) - SiabcBBBB(i,c,b,a);
    }
  }

}

template<>
void MOIntegrals<double>::formIJKA(bool doDBar){
  if(this->haveMOijka && (this->ijkaIsDBar == doDBar)) return;
  else if(this->ijkaIsDBar != doDBar) {
    if(this->Ref_ == SingleSlater<double>::TCS)
      this->ijka_.reset();
    else {
      this->ijkaAAAA_.reset();
      this->ijkaAABB_.reset();
      if(!this->singleSlater_->isClosedShell) 
        this->ijkaBBBB_.reset();
    }
  }

  this->getLocMO();

  if(this->Ref_ == SingleSlater<double>::TCS){
    this->ijka_ = std::unique_ptr<RealTensor4d>(
      new RealTensor4d(this->nO_,this->nO_,this->nO_,this->nV_));
  } else {
    this->ijkaAAAA_ = std::unique_ptr<RealTensor4d>(
      new RealTensor4d(this->nOA_,this->nOA_,this->nOA_,this->nVA_));
    this->ijkaAABB_ = std::unique_ptr<RealTensor4d>(
      new RealTensor4d(this->nOA_,this->nOA_,this->nOB_,this->nVB_));
    if(!this->singleSlater_->isClosedShell)
      this->ijkaBBBB_ = std::unique_ptr<RealTensor4d>(
        new RealTensor4d(this->nOB_,this->nOB_,this->nOB_,this->nVB_));
  }

  RealTensor4d Iinls , Iijls  , Iijks   ;
  RealTensor4d Sijka ;
  RealTensor4d IinlsA, IijlsAA, IijksAAA;  
  RealTensor4d IinlsB, IijlsBB, IijksBBB;
  RealTensor4d IijksAAB, SijkaAAAA, SijkaBBBB;
  auto NTCSxNBASIS = this->nTCS_*this->nBasis_;

  if(this->Ref_ == SingleSlater<double>::TCS){
    Iinls = RealTensor4d(this->nO_,NTCSxNBASIS,NTCSxNBASIS,NTCSxNBASIS);
    Iijls = RealTensor4d(this->nO_,this->nO_,NTCSxNBASIS,NTCSxNBASIS);
    Iijks = RealTensor4d(this->nO_,this->nO_,this->nO_,NTCSxNBASIS);
    if(doDBar)
      Sijka  = RealTensor4d(this->nO_,this->nO_,this->nO_,this->nV_);
  } else {
    IinlsA    = RealTensor4d(this->nOA_,this->nBasis_,this->nBasis_,this->nBasis_);
    IijlsAA   = RealTensor4d(this->nOA_,this->nOA_,this->nBasis_,this->nBasis_);
    IijksAAA  = RealTensor4d(this->nOA_,this->nOA_,this->nOA_,this->nBasis_);
    if(!this->singleSlater_->isClosedShell){
      IinlsB     = RealTensor4d(this->nOB_,this->nBasis_,this->nBasis_,this->nBasis_);
      IijlsBB    = RealTensor4d(this->nOB_,this->nOB_,this->nBasis_,this->nBasis_);
      IijksBBB   = RealTensor4d(this->nOB_,this->nOB_,this->nOB_,this->nBasis_);
      IijksAAB   = RealTensor4d(this->nOA_,this->nOA_,this->nOB_,this->nBasis_);
      if(doDBar){
        SijkaAAAA  = RealTensor4d(this->nOA_,this->nOA_,this->nOA_,this->nVA_);
        SijkaBBBB  = RealTensor4d(this->nOB_,this->nOB_,this->nOB_,this->nVB_);
      }
    }
  }

  enum{mu,nu,lm,sg,i,j,k,a};

  if(this->Ref_ == SingleSlater<double>::TCS){

    // First Quarter Transformation  (a ν | λ κ)
    contract(1.0,(*this->locMOOcc_),{mu,i},(*this->aointegrals_->aoERI_),{mu,nu,lm,sg},
             0.0,Iinls,{i,nu,lm,sg});
 
    // First Half Transformation     (a b | λ κ)
    contract(1.0,(*this->locMOOcc_),{nu,j},Iinls,{i,nu,lm,sg},0.0,Iijls,{i,j,lm,sg});
 
    // Third Quarter Transformation  (a b | c κ)
    contract(1.0,(*this->locMOOcc_),{lm,k},Iijls,{i,j,lm,sg},0.0,Iijks,{i,j,k,sg});

  } else {

    // First Quarter Transformation Alpha             (a ν | λ κ) [A]
    contract(1.0,(*this->locMOAOcc_),{mu,i},(*this->aointegrals_->aoERI_),{mu,nu,lm,sg},
             0.0,IinlsA,{i,nu,lm,sg});
 
    // First Half Transformation Alpha-Alpha          (a b | λ κ) [AA]
    contract(1.0,(*this->locMOAOcc_),{nu,j},IinlsA,{i,nu,lm,sg},0.0,IijlsAA,{i,j,lm,sg});
 
    // Third Quarter Transformation Alpha-Alpha-Alpha (a b | c κ) [AA|A]
    contract(1.0,(*this->locMOAOcc_),{lm,k},IijlsAA,{i,j,lm,sg},0.0,IijksAAA,{i,j,k,sg});
 
    /******************************/
    /* ONLY BUILD IF CLOSED SHELL */
    /******************************/
    if(!this->singleSlater_->isClosedShell){
 
      // First Quarter Transformation Beta            (a ν | λ κ) [B]
      contract(1.0,(*this->locMOBOcc_),{mu,i},(*this->aointegrals_->aoERI_),{mu,nu,lm,sg},
               0.0,IinlsB,{i,nu,lm,sg});
 
      // First Half Transformation Beta-Beta          (a b | λ κ) [BB]
      contract(1.0,(*this->locMOBOcc_),{nu,j},IinlsB,{i,nu,lm,sg},0.0,IijlsBB,{i,j,lm,sg});
 
      // Third Quarter Transformation Beta-Beta-Beta    (a b | c κ) [BB|B]
      contract(1.0,(*this->locMOBOcc_),{lm,k},IijlsBB,{i,j,lm,sg},0.0,IijksBBB,{i,j,k,sg});
 
      // Third Quarter Transformation Alpha-Alpha-Beta  (a b | c κ) [AA|B]
      contract(1.0,(*this->locMOBOcc_),{lm,k},IijlsAA,{i,j,lm,sg},0.0,IijksAAB,{i,j,k,sg});
 
    }
  }

  /******************************/
  /*  IF ONLY DOING SINGLE BAR  */
  /******************************/
  if(!doDBar){
    /*
     * If only doing single bar integrals, populate the pure spin
     * storage with the single bar MO integrals
     */ 
 
    if(this->Ref_ == SingleSlater<double>::TCS){
      // Last Quarter Transformation (a b | c d) 
      contract(1.0,(*this->locMOVir_),{sg,a},Iijks,{i,j,k,sg},
               0.0,(*this->ijka_),{i,j,k,a});
    } else {

      // Last Quarter Transformation Alpha-Alpha-Alpha-Alpha (a b | c d) [AA|AA]
      contract(1.0,(*this->locMOAVir_),{sg,a},IijksAAA,{i,j,k,sg},
               0.0,(*this->ijkaAAAA_),{i,j,k,a});
     
      /******************************/
      /* ONLY BUILD IF CLOSED SHELL */
      /******************************/
      // Last Quarter Transformation Beta-Beta-Beta-Beta (a b | c d) [BB|BB]
      if(!this->singleSlater_->isClosedShell)
        contract(1.0,(*this->locMOBVir_),{sg,a},IijksBBB,{i,j,k,sg},
                 0.0,(*this->ijkaBBBB_),{i,j,k,a});
    }

  /****************************************/
  /*  IF DOING DOUBLE AND OPEN-SHELL BAR  */
  /****************************************/
  } else if(!this->singleSlater_->isClosedShell || this->Ref_ == SingleSlater<double>::TCS) {
    /*
     * If doing double bar integrals with an open-shell or two-component reference, 
     * populate previously allocated intermetiates with the single bar MO integrals
     */ 

    if(this->Ref_ == SingleSlater<double>::TCS)
      // Last Quarter Transformation (a b | c d)
      contract(1.0,(*this->locMOVir_),{sg,a},Iijks,{i,j,k,sg},0.0,Sijka,{i,j,k,a});
    else {
      // Last Quarter Transformation Alpha-Alpha-Alpha-Alpha (a b | c d) [AA|AA]
      contract(1.0,(*this->locMOAVir_),{sg,a},IijksAAA,{i,j,k,sg},0.0,SijkaAAAA,{i,j,k,a});
      // Last Quarter Transformation Beta-Beta-Beta-Beta (a b | c d) [BB|BB]
      contract(1.0,(*this->locMOBVir_),{sg,a},IijksBBB,{i,j,k,sg},0.0,SijkaBBBB,{i,j,k,a});
    }
  }

  if(this->Ref_ != SingleSlater<double>::TCS){
    /*
     * Regardless of whether or not we are doing double-bar integrals, populate
     * the mixed spin storage with single-bar integrals (for RHF and UHF, the spin
     * integration removes the exchange term)
     */ 
    if(!this->singleSlater_->isClosedShell){
      // (UHF) Last Quarter Transformation Alpha-Alpha-Beta-Beta (a b | c d) [AA|BB]
      contract(1.0,(*this->locMOBVir_),{sg,a},IijksAAB,{i,j,k,sg},
               0.0,(*this->ijkaAABB_),{i,j,k,a});
    } else {
      // (RHF) Last Quarter Transformation Alpha-Alpha-Beta-Beta (a b | c d) [AA|BB]
      contract(1.0,(*this->locMOAVir_),{sg,a},IijksAAA,{i,j,k,sg},
               0.0,(*this->ijkaAABB_),{i,j,k,a});
    }
  }
  if(doDBar){
    /// THIS IS WRONG --- DO NOT USE
    if(this->Ref_ == SingleSlater<double>::TCS)
      for(auto i = 0; i < this->nO_; i++)
      for(auto j = 0; j < this->nO_; j++)
      for(auto k = 0; k < this->nO_; k++)
      for(auto a = 0; a < this->nV_; a++)
        (*this->ijka_)(i,j,k,a) = Sijka(i,j,k,a) - Sijka(i,a,k,j);
    else {
      for(auto i = 0; i < this->nOA_; i++)
      for(auto j = 0; j < this->nOA_; j++)
      for(auto k = 0; k < this->nOA_; k++)
      for(auto a = 0; a < this->nVA_; a++){
        /*
         * In the case of RHF, we can reuse the "mixed-spin" storage to build the
         * double bar integrals, as the mixed-spin storage is simply the single
         * bar integrals. This does not work for UHF
         */ 
        if(this->singleSlater_->isClosedShell)
          (*this->ijkaAAAA_)(i,j,k,a) =
            (*this->ijkaAABB_)(i,j,k,a)-(*this->ijkaAABB_)(i,a,k,j);
        else {
          (*this->ijkaAAAA_)(i,j,k,a) = SijkaAAAA(i,j,k,a) - SijkaAAAA(i,a,k,j);
          (*this->ijkaBBBB_)(i,j,k,a) = SijkaBBBB(i,j,k,a) - SijkaBBBB(i,a,k,j);
        }
      }
      if(!this->singleSlater_->isClosedShell)
        for(auto i = 0; i < this->nOB_; i++)
        for(auto j = 0; j < this->nOB_; j++)
        for(auto k = 0; k < this->nOB_; k++)
        for(auto a = 0; a < this->nVB_; a++)
          (*this->ijkaBBBB_)(i,j,k,a) = SijkaBBBB(i,j,k,a) - SijkaBBBB(i,a,k,j);
    }
  }
}

template<>
void MOIntegrals<double>::formDBar(){

  this->formIJAB(false);
  this->formIJKL(false);
  this->formABCD(false);
  this->formIAJB(false);
  this->formIABC(false);
  this->formIJKA(false);

  RealTensor4d Sijka, Sijab, Siabc, Siajb, Sabcd, Sijkl;
  RealTensor4d Dijka, Dijab, Diabc, Diajb, Dabcd, Dijkl;
  Sijka = RealTensor4d(*this->ijka_);
  Sijab = RealTensor4d(*this->ijab_);
  Siabc = RealTensor4d(*this->iabc_);
  Siajb = RealTensor4d(*this->iajb_);
  Sabcd = RealTensor4d(*this->abcd_);
  Sijkl = RealTensor4d(*this->ijkl_);

  if(this->Ref_ == SingleSlater<double>::TCS) {
    double SUM = 0.0;

    // IJKA
    for(auto i = 0; i < this->nO_; i++)
    for(auto j = 0; j < this->nO_; j++)
    for(auto k = 0; k < this->nO_; k++)
    for(auto a = 0; a < this->nV_; a++) {
      (*this->ijka_)(i,j,k,a) = Sijka(i,k,j,a) - Sijka(j,k,i,a);
    }

    Dijka = RealTensor4d(*this->ijka_);
    SUM = 0.0;
    for(double x : Dijka) SUM += x*x;
    cout << std::setprecision(12) << "IJKA = " << SUM << endl;
    

    // IJAB
    for(auto i = 0; i < this->nO_; i++)
    for(auto j = 0; j < this->nO_; j++)
    for(auto a = 0; a < this->nV_; a++)
    for(auto b = 0; b < this->nV_; b++) {
      (*this->ijab_)(i,j,a,b) = Siajb(i,a,j,b) - Siajb(i,b,j,a);
    }

    Dijab = RealTensor4d(*this->ijab_);
    SUM = 0.0;
    for(double x : Dijab) SUM += x*x;
    cout << "IJAB = " << SUM << endl;

    // ABCD
    for(auto a = 0; a < this->nV_; a++)
    for(auto b = 0; b < this->nV_; b++)
    for(auto c = 0; c < this->nV_; c++)
    for(auto d = 0; d < this->nV_; d++) {
      (*this->abcd_)(a,b,c,d) = Sabcd(a,c,b,d) - Sabcd(a,d,b,c);
    }

    Dabcd = RealTensor4d(*this->abcd_);
    SUM = 0.0;
    for(double x : Dabcd) SUM += x*x;
    cout << "ABCD = " << SUM << endl;

    // IAJB
    for(auto i = 0; i < this->nO_; i++)
    for(auto a = 0; a < this->nV_; a++)
    for(auto j = 0; j < this->nO_; j++)
    for(auto b = 0; b < this->nV_; b++) {
      (*this->iajb_)(i,a,j,b) = Sijab(i,j,a,b) - Siajb(i,b,j,a);
    }
  
    Diajb = RealTensor4d(*this->iajb_);
    SUM = 0.0;
    for(double x : Diajb) SUM += x*x;
    cout << "IAJB = " << SUM << endl;

    // IABC
    for(auto i = 0; i < this->nO_; i++)
    for(auto a = 0; a < this->nV_; a++)
    for(auto b = 0; b < this->nV_; b++)
    for(auto c = 0; c < this->nV_; c++) {
      (*this->iabc_)(i,a,b,c) = Siabc(i,b,a,c) - Siabc(i,c,a,b);
    }

    Diabc = RealTensor4d(*this->iabc_);
    SUM = 0.0;
    for(double x : Diabc) SUM += x*x;
    cout << "IABC = " << SUM << endl;

    // IJKL
    for(auto i = 0; i < this->nO_; i++)
    for(auto j = 0; j < this->nO_; j++)
    for(auto k = 0; k < this->nO_; k++)
    for(auto l = 0; l < this->nO_; l++) {
      (*this->ijkl_)(i,j,k,l) = Sijkl(i,k,j,l) - Sijkl(i,l,j,k);
    }

    Dijkl = RealTensor4d(*this->ijkl_);
    SUM = 0.0;
    for(double x : Dijkl) SUM += x*x;
    cout << "IJKL = " << SUM << endl;
  }
}

}//namespace ChronusQ
