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
#include <quantum.h>
namespace ChronusQ {
  template<>
  template<>
  Quantum<double>::Quantum(const Quantum<dcomplex> &other){
    this->elecDipole_            = 
      const_cast<Quantum<dcomplex>&>(other).elecDipole();
    this->elecQuadpole_          = 
      const_cast<Quantum<dcomplex>&>(other).elecQuadpole();
    this->elecTracelessQuadpole_ = 
      const_cast<Quantum<dcomplex>&>(other).elecTracelessQuadpole();
    this->elecOctpole_           = 
      const_cast<Quantum<dcomplex>&>(other).elecOctpole();
    this->nTCS_                  = 
      const_cast<Quantum<dcomplex>&>(other).nTCS();
    this->isClosedShell          = 
      const_cast<Quantum<dcomplex>&>(other).isClosedShell;
    this->maxMultipole_          = 
      const_cast<Quantum<dcomplex>&>(other).maxMultipole(); 

    this->isAllocated_ = false;
    auto NB  = const_cast<Quantum<dcomplex>&>(other).onePDMScalar()->rows(); 
    this->alloc(NB);

    for(auto IDen = 0; IDen < this->onePDM().size(); IDen++)
      *this->onePDM_[IDen] = 
        const_cast<Quantum<dcomplex>&>(other).onePDM()[IDen]->real();

/*
    this->isScattered_ = 
      const_cast<Quantum<dcomplex>&>(other).isScattered();
*/

/*
    auto NBT = const_cast<Quantum<dcomplex>&>(other).onePDMA()->rows(); 
    auto NB  = const_cast<Quantum<dcomplex>&>(other).onePDMScalar()->rows(); 
    auto NBTSq = NBT*NBT;
    auto NBSq = NB*NB; 
*/

/*
    this->onePDMA_ = std::unique_ptr<RealMap>(
        new RealMap( this->memManager_->malloc<double>(NBTSq),NBT,NBT)
      );
    (*this->onePDMA_) = 
      const_cast<Quantum<dcomplex>&>(other).onePDMA()->real();
    
    if(!this->isClosedShell && this->nTCS_ != 2){
      this->onePDMB_ = std::unique_ptr<RealMap>(
          new RealMap( this->memManager_->malloc<double>(NBTSq),NBT,NBT)
        );
      (*this->onePDMB_) = 
        const_cast<Quantum<dcomplex>&>(other).onePDMB()->real();
    }

    this->onePDMScalar_ = std::unique_ptr<RealMap>(
        new RealMap( this->memManager_->malloc<double>(NBSq),NB,NB)
      );
    (*this->onePDMScalar_) = 
      const_cast<Quantum<dcomplex>&>(other).onePDMScalar()->real();

    if(!this->isClosedShell || this->nTCS_ == 2) {
      this->onePDMMz_ = std::unique_ptr<RealMap>(
          new RealMap( this->memManager_->malloc<double>(NBSq),NB,NB)
        );
      (*this->onePDMMz_) = 
        const_cast<Quantum<dcomplex>&>(other).onePDMMz()->real();
    }

    if(this->nTCS_ == 2) {
      this->onePDMMx_ = std::unique_ptr<RealMap>(
          new RealMap( this->memManager_->malloc<double>(NBSq),NB,NB)
        );

      this->onePDMMy_ = std::unique_ptr<RealMap>(
          new RealMap( this->memManager_->malloc<double>(NBSq),NB,NB)
        );


      (*this->onePDMMx_) = 
        const_cast<Quantum<dcomplex>&>(other).onePDMMx()->real();
      (*this->onePDMMy_) = 
        const_cast<Quantum<dcomplex>&>(other).onePDMMy()->real();
    };
*/
  };

  template<>
  template<>
  Quantum<dcomplex>::Quantum(const Quantum<double> &other){
    this->elecDipole_            = 
      const_cast<Quantum<double>&>(other).elecDipole();
    this->elecQuadpole_          = 
      const_cast<Quantum<double>&>(other).elecQuadpole();
    this->elecTracelessQuadpole_ = 
      const_cast<Quantum<double>&>(other).elecTracelessQuadpole();
    this->elecOctpole_           = 
      const_cast<Quantum<double>&>(other).elecOctpole();
    this->nTCS_                  = 
      const_cast<Quantum<double>&>(other).nTCS();
    this->isClosedShell          = 
      const_cast<Quantum<double>&>(other).isClosedShell;
    this->maxMultipole_          = 
      const_cast<Quantum<double>&>(other).maxMultipole(); 
    this->memManager_ =
      const_cast<Quantum<double>&>(other).memManager();

    this->isAllocated_ = false;
    auto NB  = const_cast<Quantum<double>&>(other).onePDMScalar()->rows(); 
    this->alloc(NB);

    for(auto IDen = 0; IDen < this->onePDM().size(); IDen++){
      this->onePDM_[IDen]->real() = 
        *const_cast<Quantum<double>&>(other).onePDM()[IDen];
    }

/*
    this->isScattered_ = 
      const_cast<Quantum<double>&>(other).isScattered();


    auto NBT = const_cast<Quantum<double>&>(other).onePDMA()->rows(); 
    this->alloc(NBT/this->nTCS_);

    this->onePDMA_->real() = (*const_cast<Quantum<double>&>(other).onePDMA());

    if(!this->isClosedShell || this->nTCS_ == 2){
      this->onePDMScalar_->real() = 
        (*const_cast<Quantum<double>&>(other).onePDMScalar());
      this->onePDMMz_->real() = 
        (*const_cast<Quantum<double>&>(other).onePDMMz());
      if(this->nTCS_ == 1)
        this->onePDMB_->real() = 
          (*const_cast<Quantum<double>&>(other).onePDMB());
      else {
        this->onePDMMx_->real() = 
          (*const_cast<Quantum<double>&>(other).onePDMMx());
        this->onePDMMy_->real() = 
          (*const_cast<Quantum<double>&>(other).onePDMMy());

      } // NTCS check
    } // if not RHF/KS
*/
  }; // COPY
};
