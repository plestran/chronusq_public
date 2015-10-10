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
  bool doTCS = (this->Ref_ == TCS);
  bool doRHF = (this->isClosedShell && !doTCS);
  if(!this->haveDensity) this->formDensity();
  if(this->controls_->directTwoE && !this->controls_->doDF)
    this->aointegrals_->twoEContractDirect(doRHF,true,false,doTCS,*this->densityA_,*this->PTA_,*this->densityB_,*this->PTB_);
  else if(this->controls_->doDF)
    this->aointegrals_->twoEContractDF(doRHF,true,*this->densityA_,*this->PTA_,*this->densityB_,*this->PTB_);
  else
    this->aointegrals_->twoEContractN4(doRHF,true,false,doTCS,*this->densityA_,*this->PTA_,*this->densityB_,*this->PTB_);
  if(this->controls_->printLevel >= 3) this->printPT();
//if(doTCS)CErr();
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
//if(this->Ref_ == TCS) this->basisset_->resetMapSh2Bf();
  if(!this->aointegrals_->haveAOOneE) this->aointegrals_->computeAOOneE();
//if(this->Ref_ == TCS) {
//  this->basisset_->resetMapSh2Bf(); 
//  this->basisset_->makeMapSh2Bf(this->nTCS_);
//}
// Form Vxc for DFT
  if(this->controls_->DFT) this->formVXC();
  this->fockA_->setZero();
/*
  if(this->Ref_ != TCS) fockA_->real()+=(*this->aointegrals_->oneE_);
  else {
    for(auto I = 0, i = 0; i < this->nBasis_; I += 2, i++)
    for(auto J = 0, j = 0; j < this->nBasis_; J += 2, j++){
      this->fockA_->real()(I,J)     += (*this->aointegrals_->oneE_)(i,j);
      this->fockA_->real()(I+1,J+1) += (*this->aointegrals_->oneE_)(i,j);
    }
  }
*/
  fockA_->real()+=(*this->aointegrals_->oneE_);
#ifndef USE_LIBINT
  fockA_->real()+=(*this->coulombA_);
  fockA_->real()-=(*this->exchangeA_);
#else
  *(fockA_)+=(*this->PTA_);
#endif
  if(this->controls_->DFT){
    cout << "Single Slater Numeric : Print" <<endl;
    cout << (*this->vXCA_)  << endl;
    (*this->fockA_) += (*this->vXCA_);
  }
  if(!this->isClosedShell && this->Ref_ != TCS){
    this->fockB_->setZero();
    fockB_->real()+=(*this->aointegrals_->oneE_);
#ifndef USE_LIBINT
    fockB_->real()+=(*this->coulombB_);
    fockB_->real()-=(*this->exchangeB_);
#else
    *(fockB_)+=(*this->PTB_);
#endif
    if(this->controls_->DFT){
    cout << "Single Slater Numeric : Print" <<endl;
    cout << (*this->vXCB_)  << endl;
      (*this->fockB_) += (*this->vXCB_);
    }
  };

  // Add in the electric field component if they are non-zero
  std::array<double,3> null{{0,0,0}};
  if(this->elecField_ != null){
    //this->fileio_->out << "Adding in Electric Field Contribution" << endl;
    int NB = this->nTCS_*this->nBasis_;
    int NBSq = NB*NB;
    int iBuf = 0;
    for(auto iXYZ = 0; iXYZ < 3; iXYZ++){
      ConstRealMap mu(&this->aointegrals_->elecDipole_->storage()[iBuf],NB,NB);
      fockA_->real() += this->elecField_[iXYZ] * mu;
      if(!this->isClosedShell && this->Ref_ != TCS) 
        fockB_->real() += this->elecField_[iXYZ] * mu;
      iBuf += NBSq;
    }
  }
  if(this->controls_->printLevel>=2) this->printFock(); 
};

