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
// Only works for RHF FIXME
// Only available for Libint Integrals
#ifdef USE_LIBINT
template<typename T>
void SingleSlater<T>::formPT(){
  if(!this->haveDensity) this->formDensity();
  if(this->controls_->directTwoE)
    this->aointegrals_->twoEContractDirect(true,*this->densityA_,*this->PTA_);
  else if(this->controls_->doDF)
    this->aointegrals_->twoEContractDF(true,*this->densityA_,*this->PTA_);
  else
    this->aointegrals_->twoEContractN4(true,*this->densityA_,*this->PTA_);
  if(this->controls_->printLevel >= 3) prettyPrint(this->fileio_->out,(*this->PTA_),"Alpha Perturbation Tensor");
}
#endif

/********************
 * Form Fock Matrix *
 ********************/
// Only reliable for RHF (b/c factors of 2) FIXME
// In-house integrals are broken for Fock   FIXME
template<typename T>
void SingleSlater<T>::formFock(){
  if(!this->haveDensity) this->formDensity();
#ifndef USE_LIBINT
  if(!this->haveCoulomb) this->formCoulomb();
  if(!this->haveExchange) this->formExchange();
#else
  if(!this->havePT) this->formPT();
#endif
  if(!this->aointegrals_->haveAOOneE) this->aointegrals_->computeAOOneE();

  this->fockA_->setZero();
  *(fockA_)+=(*this->aointegrals_->oneE_);
#ifndef USE_LIBINT
  *(fockA_)+=2*(*this->coulombA_);
  *(fockA_)-=2*(*this->exchangeA_);
#else
  *(fockA_)+=2*(*this->PTA_);
#endif
  if(!this->RHF_){
    this->fockB_->setZero();
    *(fockB_)+=(*this->aointegrals_->oneE_);
#ifndef USE_LIBINT
    *(fockB_)+=2*(*this->coulombB_);
    *(fockB_)-=2*(*this->exchangeB_);
#else
    *(fockB_)+=(*this->PTB_);
#endif
  };
  if(this->controls_->printLevel>=2) {
    prettyPrint(this->fileio_->out,(*this->fockA_),"Alpha Fock");
    if(!this->RHF_) prettyPrint(this->fileio_->out,(*this->fockB_),"Beta Fock");
  };
};

