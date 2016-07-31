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
  else if(this->guess_ == RANDOM) this->RandomGuess();
  else CErr("Guess NYI",this->fileio_->out);
}
template<typename T>
void SingleSlater<T>::COREGuess(){
  this->aointegrals_->computeAOOneE(); 
  if(this->nTCS_ == 1 && this->isClosedShell){
    this->fockA_->real() = (*this->aointegrals_->coreH_);
  } else {
    this->fockScalar_->real() = 2*(*this->aointegrals_->coreH_);
    this->fockMz_->setZero();

    std::vector<std::reference_wrapper<TMap>> toGather;
    toGather.emplace_back(*this->fockScalar_);
    toGather.emplace_back(*this->fockMz_);

    if(this->nTCS_ ==  2) {
      this->fockMy_->setZero();
      this->fockMx_->setZero();
      toGather.emplace_back(*this->fockMy_);
      toGather.emplace_back(*this->fockMx_);

      Quantum<T>::spinGather(*this->fockA_,toGather);
    } else 
      Quantum<T>::spinGather(*this->fockA_,*this->fockB_,toGather);
  }
  this->haveMO = true;
  this->haveDensity = true;

  //this->orthoFock();
  this->orthoFock3();
  if(this->printLevel_ > 3){
    prettyPrint(this->fileio_->out,*this->fockA_,"Initial FA");
    prettyPrint(this->fileio_->out,*this->fockOrthoA_,"Initial FOrthoA");
  }
  this->diagFock2();
  this->formDensity();
  this->cpyAOtoOrthoDen();
  this->unOrthoDen();
  if(this->printLevel_ > 3){
    prettyPrint(this->fileio_->out,*this->onePDMA_,"Initial PA");
    prettyPrint(this->fileio_->out,*this->onePDMOrthoA_,"Initial POrthoA");
  }
  this->computeEnergy();
  this->backTransformMOs();
  if(this->printLevel_ > 3) {
    prettyPrint(this->fileio_->out,*this->moA_,"Initial MOA");
  }

};
