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
template<typename T>
void SingleSlater<T>::formGuess(){
  if(this->guess_ == SAD) this->SADGuess();
  else if(this->guess_ == CORE) this->COREGuess();
  else if(this->guess_ == READ) this->READGuess();
  else CErr("Guess NYI",this->fileio_->out);
}
template<typename T>
void SingleSlater<T>::COREGuess(){
 
  this->onePDMA_->setZero();
  this->moA_->setZero();
  if(this->nTCS_ == 2 || !this->isClosedShell){
    if(this->nTCS_ == 1) this->onePDMB_->setZero();
    this->onePDMScalar_->setZero();
    this->onePDMMz_->setZero();
    if(this->nTCS_ == 2) {
      this->onePDMMy_->setZero();
      this->onePDMMx_->setZero();
    }
  }
  this->haveMO = true;
  this->haveDensity = true;
};
