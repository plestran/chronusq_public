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
  bool doTCS = (this->Ref_ == TCS);
  bool doRHF = (this->isClosedShell && !doTCS);
  bool doKS  = this->isDFT;
  if(!this->haveDensity) this->formDensity();


  if(this->aointegrals_->integralAlgorithm == AOIntegrals::DIRECT)
  {; }
//  this->aointegrals_->twoEContractDirect(doRHF,doKS,true,false,doTCS,
//  *this->onePDMA_,*this->PTA_,*this->onePDMB_,*this->PTB_);
  else if(this->aointegrals_->integralAlgorithm == AOIntegrals::DENFIT)
    this->aointegrals_->twoEContractDF(doRHF,doKS,true,*this->onePDMA_,
    *this->PTA_,*this->onePDMB_,*this->PTB_);
  else if(this->aointegrals_->integralAlgorithm == AOIntegrals::INCORE)
    this->aointegrals_->twoEContractN4(doRHF,doKS,true,false,doTCS,
    *this->onePDMA_,*this->PTA_,*this->onePDMB_,*this->PTB_);
  if(this->printLevel_ >= 3 && getRank() == 0) this->printPT();

  std::vector<std::reference_wrapper<TMatrix>> mats;
  std::vector<std::reference_wrapper<TMatrix>> ax;
  std::vector<AOIntegrals::ERI_CONTRACTION_TYPE> contList;
  std::vector<double> scalingFactors;
  this->scatterDensity();

  if(this->nTCS_ == 1 && this->isClosedShell) {
    TMatrix GPScalar(this->nBasis_,this->nBasis_);

    mats.emplace_back(*this->onePDMA_);
    mats.emplace_back(*this->onePDMA_);
    ax.emplace_back(GPScalar);
    ax.emplace_back(GPScalar);
    contList.push_back(AOIntegrals::ERI_CONTRACTION_TYPE::COULOMB);
    contList.push_back(AOIntegrals::ERI_CONTRACTION_TYPE::EXCHANGE);

    scalingFactors.push_back(1.0);
    scalingFactors.push_back(-0.5);

    if(this->aointegrals_->integralAlgorithm == AOIntegrals::DIRECT){
      this->aointegrals_->newTwoEContractDirect(mats,ax,contList,
          scalingFactors);
    
      (*this->PTA_) = GPScalar * 0.5;
    }

  } else {
    TMatrix GPScalar(this->nBasis_,this->nBasis_);
    TMatrix GPMz(this->nBasis_,this->nBasis_);
   
    mats.emplace_back(*this->onePDMScalar_);
    mats.emplace_back(*this->onePDMScalar_);
    mats.emplace_back(*this->onePDMMz_);
    ax.emplace_back(GPScalar);
    ax.emplace_back(GPScalar);
    ax.emplace_back(GPMz);
    contList.push_back(AOIntegrals::ERI_CONTRACTION_TYPE::COULOMB);
    contList.push_back(AOIntegrals::ERI_CONTRACTION_TYPE::EXCHANGE);
    contList.push_back(AOIntegrals::ERI_CONTRACTION_TYPE::EXCHANGE);
    scalingFactors.push_back(1.0);
    scalingFactors.push_back(-0.5);
    scalingFactors.push_back(-0.5);

    if(this->aointegrals_->integralAlgorithm == AOIntegrals::DIRECT){
      this->aointegrals_->newTwoEContractDirect(mats,ax,contList,
          scalingFactors);
    
      std::vector<std::reference_wrapper<TMatrix>> toGather;
      toGather.push_back(GPScalar);
      toGather.push_back(GPMz);
      Quantum<T>::spinGather((*this->PTA_),(*this->PTB_),toGather);
    }
  }
  this->gatherDensity();
}
#endif

/********************
 * Form Fock Matrix *
 ********************/
template<typename T>
void SingleSlater<T>::formFock(){
#ifdef CQ_ENABLE_MPI
  // Syncronize MPI processes before fock build
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  
  if(!this->haveDensity) this->formDensity();
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

    if (this->isDFT){
//    Timing
      std::chrono::high_resolution_clock::time_point start;
      std::chrono::high_resolution_clock::time_point finish;
      std::chrono::duration<double> duration_formVxc;
     
      if(this->printLevel_ >= 3) {
        start = std::chrono::high_resolution_clock::now();
      }

//    this->formVXC();
      this->formVXC_store();
     
      if(this->printLevel_ >= 3) {
        finish = std::chrono::high_resolution_clock::now();
        duration_formVxc = finish - start;
        this->fileio_->out << endl << "CPU time for VXC integral:  "
                           << duration_formVxc.count() << " seconds." 
                           << endl;
      }
    }
  
    this->fockA_->setZero();
    this->fockA_->real() += (*this->aointegrals_->oneE_);
#ifndef USE_LIBINT
    this->fockA_->real() += (*this->coulombA_);
    this->fockA_->real() -= (*this->exchangeA_);
#else
    (*this->fockA_) += (*this->PTA_);
#endif
    if(this->isDFT){ 
      (*this->fockA_) += (*this->vXA_);
      (*this->fockA_) += (*this->vCorA_);
    }

    if(!this->isClosedShell && this->Ref_ != TCS){
      this->fockB_->setZero();
      fockB_->real()+=(*this->aointegrals_->oneE_);
#ifndef USE_LIBINT
      fockB_->real()+=(*this->coulombB_);
      fockB_->real()-=(*this->exchangeB_);
#else
      *(fockB_) += (*this->PTB_);
#endif
      if(this->isDFT) {
        (*fockB_) += (*this->vXB_);
        (*fockB_) += (*this->vCorB_);
      }
    }
    // Add in the electric field component if they are non-zero
    std::array<double,3> null{{0,0,0}};
    if(this->elecField_ != null){
      //this->fileio_->out << "Adding in Electric Field Contribution" << endl;
      int NB = this->nTCS_*this->nBasis_;
      int NBSq = NB*NB;
      int iBuf = 0;
      for(auto iXYZ = 0; iXYZ < 3; iXYZ++){
        ConstRealMap mu(&this->aointegrals_->elecDipole_->storage()[iBuf],
          NB,NB);
        fockA_->real() += this->elecField_[iXYZ] * mu;
        if(!this->isClosedShell && this->Ref_ != TCS) 
          fockB_->real() += this->elecField_[iXYZ] * mu;
        iBuf += NBSq;
      }
    }
    if(this->printLevel_ >= 2) this->printFock(); 
  }
#ifdef CQ_ENABLE_MPI
  // Syncronize MPI processes after fock build
  MPI_Barrier(MPI_COMM_WORLD);
#endif
};

