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
template <typename T>
void WaveFunction<T>::alloc() {
  this->checkMeta();
  Quantum<T>::alloc(this->nBasis_);

  auto NB = this->nBasis_;
  auto NBT = this->nTCS_ * NB;
  auto NBSq = NB*NB;
  auto NBTSq = NBT*NBT;

  // MO Coeffs
  this->moA_ = std::unique_ptr<TMap>(
    new TMap(this->memManager_->template malloc<T>(NBTSq),NBT,NBT)); 

  this->moA_->setZero();
  if(this->nTCS_ == 1 and !this->isClosedShell) {
    this->moB_ = std::unique_ptr<TMap>(
      new TMap(this->memManager_->template malloc<T>(NBTSq),NBT,NBT)); 

    this->moB_->setZero();
  }

  // MO Eigenenergies
  this->epsA_ = std::unique_ptr<RealMap>(
    new RealMap(this->memManager_->template malloc<double>(NBTSq),NBT,1)); 
    this->epsA_->setZero();
  if(this->nTCS_ == 1 and !this->isClosedShell){
    this->epsB_ = std::unique_ptr<RealMap>(
      new RealMap(this->memManager_->template malloc<double>(NBTSq),NBT,1)); 
    this->epsB_->setZero();
  }
}
