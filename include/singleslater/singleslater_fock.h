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

/********************************
 * Form Perturbation Tensor (G) *
 ********************************/
// Only available for Libint Integrals
#ifdef USE_LIBINT
template<typename T>
void SingleSlater<T>::formPT(){

  std::vector<std::reference_wrapper<TMap>> mats;
  std::vector<std::reference_wrapper<TMap>> ax;
  std::vector<AOIntegrals::ERI_CONTRACTION_TYPE> contList;
  std::vector<double> scalingFactors;
  double exchFactor = -0.5;
  if(this->isDFT) exchFactor = 0.0;

  mats.emplace_back(*this->onePDMScalar_);
  mats.emplace_back(*this->onePDMScalar_);
  ax.emplace_back(*this->PTScalar_);
  ax.emplace_back(*this->PTScalar_);
  contList.push_back(AOIntegrals::ERI_CONTRACTION_TYPE::COULOMB);
  contList.push_back(AOIntegrals::ERI_CONTRACTION_TYPE::EXCHANGE);
  scalingFactors.push_back(1.0);
  scalingFactors.push_back(exchFactor);

  this->PTScalar_->setZero();

  if(this->nTCS_ == 2 or !this->isClosedShell) {
    mats.emplace_back(*this->onePDMMz_);
    ax.emplace_back(*this->PTMz_);
    contList.push_back(AOIntegrals::ERI_CONTRACTION_TYPE::EXCHANGE);
    scalingFactors.push_back(exchFactor);
    this->PTMz_->setZero();
  }

  if(this->nTCS_ == 2){
    mats.emplace_back(*this->onePDMMy_);
    mats.emplace_back(*this->onePDMMx_);
    ax.emplace_back(*this->PTMy_);
    ax.emplace_back(*this->PTMx_);
    contList.push_back(AOIntegrals::ERI_CONTRACTION_TYPE::EXCHANGE);
    contList.push_back(AOIntegrals::ERI_CONTRACTION_TYPE::EXCHANGE);
    scalingFactors.push_back(exchFactor);
    scalingFactors.push_back(exchFactor);
    this->PTMy_->setZero();
    this->PTMx_->setZero();
  }

  this->aointegrals_->newTwoEContract(mats,ax,contList,scalingFactors);

/*
  if(this->nTCS_ == 1 && this->isClosedShell) {
    (*this->PTA_) *= 0.5; // This is effectively a gather operation
  } else {
    std::vector<std::reference_wrapper<TMap>> toGather;
    toGather.emplace_back(*this->PTScalar_);
    toGather.emplace_back(*this->PTMz_);
    if(this->nTCS_ == 1)
      Quantum<T>::spinGather((*this->PTA_),(*this->PTB_),toGather);
    else {
      toGather.emplace_back(*this->PTMy_);
      toGather.emplace_back(*this->PTMx_);
      Quantum<T>::spinGather((*this->PTA_),toGather);
    }

  }
*/

  if(this->printLevel_ >= 3 && getRank() == 0) this->printPT();
//prettyPrint(cout,*this->onePDMA_,"P in PT");
//prettyPrint(cout,*this->PTA_,"PT in PT");
}
#endif

/********************
 * Form Fock Matrix *
 ********************/
template<typename T>
void SingleSlater<T>::formFock(bool increment){
#ifdef CQ_ENABLE_MPI
  // Syncronize MPI processes before fock build
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  
//if(!this->haveDensity) this->formDensity();
#ifndef USE_LIBINT
  if(getRank() == 0){
    // Not even sure if these guys still work let alone are MPI
    // capable
    if(!this->haveCoulomb) this->formCoulomb();
    if(!this->haveExchange) this->formExchange();
  }
#else
  // All MPI processes go to FormPT
  this->formPT();
#endif

  if(getRank() == 0) {
    if(!this->aointegrals_->haveAOOneE) this->aointegrals_->computeAOOneE();

    if (this->isDFT) this->formVXC_new();
  
/*
    if(this->nTCS_ == 1 && this->isClosedShell) {
      if(!increment) {
        cout << "Not Incrementing" << endl;
        this->fockA_->setZero();
        this->fockA_->real() += (*this->aointegrals_->coreH_);
        this->aointegrals_->addElecDipole(*this->fockA_,this->elecField_);
      }
      (*this->fockA_) += (*this->PTA_);
      if(this->isDFT){ 
        (*this->fockA_) += (*this->vXA_);
      }
    } else {
      if(!increment) {
        this->fockScalar_->setZero();
        this->fockMz_->setZero();
        this->fockScalar_->real() += (*this->aointegrals_->coreH_);
        
        if(this->nTCS_ == 2) {
            this->fockMx_->setZero();
            this->fockMy_->setZero();
          if(this->aointegrals_->doX2C){
            this->fockMx_->real() = 2*(*this->aointegrals_->oneEmx_);
            this->fockMy_->real() = 2*(*this->aointegrals_->oneEmy_);
            this->fockMz_->real() = 2*(*this->aointegrals_->oneEmz_);
            // -----------------------------------
            // SCALE SO parts by 'i'
            // -----------------------------------
            Quantum<T>::complexMyScale(*this->fockMx_);
            Quantum<T>::complexMyScale(*this->fockMy_);
            Quantum<T>::complexMyScale(*this->fockMz_);
            // -----------------------------------
          }
        }

        this->aointegrals_->addElecDipole(*this->fockScalar_,this->elecField_);
        (*this->fockScalar_) *= 2.0;
      }
      (*this->fockScalar_)      += (*this->PTScalar_);        
      (*this->fockMz_)          += (*this->PTMz_);

      // FIXME: Needs to ge generalized for the 2C case
      if(this->nTCS_ == 1 && this->isDFT){
          (*this->fockScalar_) += (*this->vXA_);
          (*this->fockScalar_) += (*this->vXB_);
          (*this->fockMz_)     += (*this->vXA_);
          (*this->fockMz_)     -= (*this->vXB_);
      }

      // Transform the Fock for CUHF
      if(this->Ref_ == CUHF) {
        this->formNO();
        this->fockCUHF();
      }

      std::vector<std::reference_wrapper<TMap>> toGather;
      toGather.emplace_back(*this->fockScalar_);
      toGather.emplace_back(*this->fockMz_);
      if(this->nTCS_ == 1)
        Quantum<T>::spinGather(*this->fockA_,*this->fockB_,toGather);
      else {
        (*this->fockMx_) += (*this->PTMx_);
        (*this->fockMy_) += (*this->PTMy_);
        toGather.emplace_back(*this->fockMy_);
        toGather.emplace_back(*this->fockMx_);
        Quantum<T>::spinGather(*this->fockA_,toGather);

      //  prettyPrint(this->fileio_->out,(*this->fockA_)," fockA_ ");
      }
    }
*/

    if(!increment) {
      this->fockScalar_->setZero();

      if(this->nTCS_ == 2 or !this->isClosedShell)
        this->fockMz_->setZero();

      if(this->nTCS_ == 2){
        this->fockMy_->setZero();
        this->fockMx_->setZero();
      }

      (*this->fockScalar_) = this->aointegrals_->coreH_->template cast<T>();
      this->aointegrals_->addElecDipole(*this->fockScalar_,this->elecField_);

      if(this->nTCS_ == 2 and this->aointegrals_->doX2C) {
        (*this->fockMx_) = this->aointegrals_->oneEmx_->template cast<T>();
        (*this->fockMy_) = this->aointegrals_->oneEmy_->template cast<T>();
        (*this->fockMz_) = this->aointegrals_->oneEmz_->template cast<T>();
        // -----------------------------------
        // SCALE SO parts by 'i'
        // -----------------------------------
        /*
        Quantum<T>::complexMyScale(*this->fockMx_);
        Quantum<T>::complexMyScale(*this->fockMy_);
        Quantum<T>::complexMyScale(*this->fockMz_);
        */
        (*this->fockMx_) *= ComplexScale<T>();
        (*this->fockMy_) *= ComplexScale<T>();
        (*this->fockMz_) *= ComplexScale<T>();
        // -----------------------------------
      }
      for(auto iF = fock_.begin(); iF != fock_.end(); iF++)
        *(*iF) *= 2;
    }

/*
      this->fileio_->out << "Before PT" << endl;
      this->printFock();
      for(auto i = 0; i < this->fockScalar_->size(); i++){
        cout << this->fockScalar_->data()[i] << " + " << this->PTScalar_->data()[i] << " = " << this->fockScalar_->data()[i] + this->PTScalar_->data()[i] << endl;
        this->fockScalar_->data()[i] += this->PTScalar_->data()[i];
      }
*/
    for(auto iF = 0; iF < fock_.size(); iF++) {
      *(fock_[iF]) += *(PT_[iF]);
      if(this->isDFT)
        (*fock_[iF]) += (*vXC_[iF]);
    }
/*
    (*this->fockScalar_) += (*this->PTScalar_);
*/
/*
      this->fileio_->out << "After PT" << endl;
      this->printFock();
*/


    if(this->printLevel_ >= 2) this->printFock(); 
  }
#ifdef CQ_ENABLE_MPI
  // Syncronize MPI processes after fock build
  MPI_Barrier(MPI_COMM_WORLD);
#endif
};

