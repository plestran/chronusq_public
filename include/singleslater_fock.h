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
/********************************
 * Form Perturbation Tensor (G) *
 ********************************/
// Only available for Libint Integrals
#ifdef USE_LIBINT
template<typename T>
void SingleSlater<T>::formPT(){
  if(!this->haveDensity) this->formDensity();
  if(this->controls_->directTwoE && !this->controls_->doDF)
    this->aointegrals_->twoEContractDirect(this->RHF_,true,false,*this->densityA_,*this->PTA_,*this->densityB_,*this->PTB_);
  else if(this->controls_->doDF)
    this->aointegrals_->twoEContractDF(this->RHF_,true,*this->densityA_,*this->PTA_,*this->densityB_,*this->PTB_);
  else
    this->aointegrals_->twoEContractN4(this->RHF_,true,*this->densityA_,*this->PTA_,*this->densityB_,*this->PTB_);
  if(this->controls_->printLevel >= 3) {
    prettyPrint(this->fileio_->out,(*this->PTA_),"Alpha Perturbation Tensor");
    if(!this->RHF_) prettyPrint(this->fileio_->out,(*this->PTB_),"Beta Perturbation Tensor");
  }
}
#endif

/********************
 * Form Fock Matrix *
 ********************/
template<typename T>
void SingleSlater<T>::formFock(){
  if(!this->haveDensity) this->formDensity();
#ifndef USE_LIBINT
  if(!this->haveCoulomb) this->formCoulomb();
  if(!this->haveExchange) this->formExchange();
#else
  this->formPT();
#endif
  if(!this->aointegrals_->haveAOOneE) this->aointegrals_->computeAOOneE();

  this->fockA_->setZero();
  fockA_->real()+=(*this->aointegrals_->oneE_);
#ifndef USE_LIBINT
  fockA_->real()+=(*this->coulombA_);
  fockA_->real()-=(*this->exchangeA_);
#else
  *(fockA_)+=(*this->PTA_);
#endif
  if(!this->RHF_){
    this->fockB_->setZero();
    fockB_->real()+=(*this->aointegrals_->oneE_);
#ifndef USE_LIBINT
    fockB_->real()+=(*this->coulombB_);
    fockB_->real()-=(*this->exchangeB_);
#else
    *(fockB_)+=(*this->PTB_);
#endif
  };
  if(this->controls_->printLevel>=2) {
    prettyPrint(this->fileio_->out,(*this->fockA_),"Alpha Fock");
    if(!this->RHF_) prettyPrint(this->fileio_->out,(*this->fockB_),"Beta Fock");
  };
};

