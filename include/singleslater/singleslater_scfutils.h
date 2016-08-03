template <typename T>
void SingleSlater<T>::initSCFPtr(){
  this->occNumMem_     = NULL;

  this->FpAlphaMem_    = NULL;
  this->FpBetaMem_     = NULL;
  this->POldAlphaMem_  = NULL;
  this->POldBetaMem_   = NULL;
  this->ErrorAlphaMem_ = NULL;
  this->ErrorBetaMem_  = NULL;
  this->FADIIS_        = NULL;
  this->FBDIIS_        = NULL;
  this->lambdaMem_     = NULL;
  this->delFMem_       = NULL;
  this->PNOMem_        = NULL;

  this->NBSQScr1_ = NULL;
  this->NBSQScr2_ = NULL;
  this->ScalarScr1_ = NULL;
  this->ScalarScr2_ = NULL;
  this->MzScr1_ = NULL;
  this->MzScr2_ = NULL;
  this->MyScr1_ = NULL;
  this->MyScr2_ = NULL;
  this->MxScr1_ = NULL;
  this->MxScr2_ = NULL;
}; // initSCFPtr

template<typename T>
void SingleSlater<T>::initSCFMem3(){
  this->initSCFPtr();

  auto NTCSxNBASIS = this->nTCS_ * this->nBasis_;
  auto NSQ = NTCSxNBASIS  * NTCSxNBASIS;

  this->POldAlphaMem_  = this->memManager_->template malloc<T>(NSQ);
  if(this->nTCS_ == 1 && !this->isClosedShell) 
    this->POldBetaMem_  = this->memManager_->template malloc<T>(NSQ);

  if(this->Ref_ == CUHF) {
    this->delFMem_   = this->memManager_->template malloc<T>(NSQ);
    this->lambdaMem_ = this->memManager_->template malloc<T>(NSQ);
    this->PNOMem_    = this->memManager_->template malloc<T>(NSQ);
    this->occNumMem_ = this->memManager_->template malloc<double>(NTCSxNBASIS);
  }

   
  if(this->doDIIS) this->initDIISFiles();
  std::vector<hsize_t> dims;
  dims.push_back(this->nBasis_);
  dims.push_back(this->nBasis_);

  this->DScalarOld_ = this->fileio_->createScratchPartition(H5PredType<T>(),
    "Most Recent Density Matrix (Scalar) for RMS",dims);

  if(this->nTCS_ == 2 || !this->isClosedShell)
    this->DMzOld_ = this->fileio_->createScratchPartition(H5PredType<T>(),
      "Most Recent Density Matrix (Mz) for RMS",dims);
  
  if(this->nTCS_ == 2) {
    this->DMyOld_ = this->fileio_->createScratchPartition(H5PredType<T>(),
      "Most Recent Density Matrix (My) for RMS",dims);
    this->DMxOld_ = this->fileio_->createScratchPartition(H5PredType<T>(),
      "Most Recent Density Matrix (Mx) for RMS",dims);
  }

};

template<typename T>
void SingleSlater<T>::cleanupSCFMem3(){

  auto NTCSxNBASIS = this->nTCS_ * this->nBasis_;
  auto NSQ = NTCSxNBASIS  * NTCSxNBASIS;
  this->memManager_->free(this->POldAlphaMem_,NSQ);  
  if(this->nTCS_ == 1 && !this->isClosedShell) 
    this->memManager_->free(this->POldBetaMem_,NSQ);  

  if(this->Ref_ == CUHF){ 
    this->memManager_->free(this->delFMem_,NSQ);   
    this->memManager_->free(this->lambdaMem_,NSQ); 
    this->memManager_->free(this->PNOMem_,NSQ);    
    this->memManager_->free(this->occNumMem_,NTCSxNBASIS); 
  }
}

template<typename T>
void SingleSlater<T>::orthoFock3(){
  if(this->nTCS_ == 1 && this->isClosedShell){
    this->aointegrals_->Ortho1Trans(
      *this->fockA_,*this->fockOrthoA_);
  } else {
    this->aointegrals_->Ortho1Trans(
      *this->fockScalar_,*this->fockOrthoScalar_);
    this->aointegrals_->Ortho1Trans(
      *this->fockMz_,*this->fockOrthoMz_);
    if(this->nTCS_ == 2){
      this->aointegrals_->Ortho1Trans(
        *this->fockMy_,*this->fockOrthoMx_);
      this->aointegrals_->Ortho1Trans(
        *this->fockMy_,*this->fockOrthoMx_);
    }
  }
  this->gatherOrthoFock();

}

template<typename T>
void SingleSlater<T>::gatherOrthoFock(){
  if(this->nTCS_ == 1 && this->isClosedShell) return;

  std::vector<std::reference_wrapper<TMap>> toGather;
  toGather.emplace_back(*this->fockOrthoScalar_);
  toGather.emplace_back(*this->fockOrthoMz_);
  if(this->nTCS_ == 1)
    Quantum<T>::spinGather(*this->fockOrthoA_,
      *this->fockOrthoB_,toGather);
  else {
    toGather.emplace_back(*this->fockOrthoMy_);
    toGather.emplace_back(*this->fockOrthoMx_);
    Quantum<T>::spinGather(*this->fockOrthoA_,toGather);
  }

};

template<typename T>
void SingleSlater<T>::unOrthoDen3(){
  if(this->nTCS_ == 1 && this->isClosedShell){
    this->aointegrals_->Ortho1Trans(
      *this->onePDMOrthoA_,*this->onePDMA_);
    prettyPrint(cout,*this->onePDMA_,"PScalar");
  } else {
    this->aointegrals_->Ortho1Trans(
      *this->onePDMOrthoScalar_,*this->onePDMScalar_);
    prettyPrint(cout,*this->onePDMScalar_,"PScalar");
    this->aointegrals_->Ortho1Trans(
      *this->onePDMOrthoMz_,*this->onePDMMz_);
    if(this->nTCS_ == 2){
      this->aointegrals_->Ortho1Trans(
        *this->onePDMOrthoMy_,*this->onePDMMx_);
      this->aointegrals_->Ortho1Trans(
        *this->onePDMOrthoMy_,*this->onePDMMx_);
    }
  }
}

template<typename T>
void SingleSlater<T>::cpyAOtoOrthoDen(){
  if(this->nTCS_ == 1 && this->isClosedShell) {
    // Just copy for RHF
    (*this->onePDMOrthoA_) = (*this->onePDMA_);
  } else if(this->nTCS_ == 2 || !this->isClosedShell){
    // Scatter 
    std::vector<std::reference_wrapper<TMap>> scattered;
    scattered.emplace_back(*this->onePDMOrthoScalar_);
    scattered.emplace_back(*this->onePDMOrthoMz_);

    if(this->nTCS_ == 1) {
      Quantum<T>::spinScatter(*this->onePDMA_,*this->onePDMB_,scattered);
    } else {
      scattered.emplace_back(*this->onePDMOrthoMy_);
      scattered.emplace_back(*this->onePDMOrthoMx_);
      Quantum<T>::spinScatter(*this->onePDMA_,scattered);
    }

    /*
    // Copy
    (*this->onePDMOrthoScalar_) = (*this->onePDMScalar_);
    (*this->onePDMOrthoMz_)     = (*this->onePDMMz_);
    if(this->nTCS_ == 2) {
      (*this->onePDMOrthoMy_)     = (*this->onePDMMy_);
      (*this->onePDMOrthoMx_)     = (*this->onePDMMx_);
    }
    */
  }
};

template<typename T>
SCFConvergence SingleSlater<T>::evalConver3(){

  if(this->nTCS_ == 1 && this->isClosedShell){
    prettyPrint(cout,*this->onePDMA_,"PScalar in Eval");
  } else {
    prettyPrint(cout,*this->onePDMScalar_,"PScalar in Eval");
  }
  // Energy Convergence
  double EOld = this->totalEnergy;
  this->formFock();
  this->computeEnergy();
  double EDelta = this->totalEnergy - EOld;

  double PARMS(0),PBRMS(0);

/*
  for(auto I = 0; I < this->onePDMA_->size(); I++){
    T DIFF = this->onePDMA_->data()[I] - this->POldAlphaMem_[I];
    DIFF = std::conj(DIFF)*DIFF;
    PARMS += reinterpret_cast<double(&)[2]>(DIFF)[0];
    if(this->nTCS_ == 1 && !this->isClosedShell) {
      DIFF = this->onePDMB_->data()[I] - this->POldBetaMem_[I];
      DIFF = std::conj(DIFF)*DIFF;
      PBRMS += reinterpret_cast<double(&)[2]>(DIFF)[0];
    }
  }
*/

  

  DScalarOld_->read(this->NBSqScratch_->data(),H5PredType<T>());
  if(this->nTCS_ == 1 && this->isClosedShell) {
    for(auto I = 0; I < this->onePDMA_->size(); I++){
      T DIFF = this->onePDMA_->data()[I] - this->NBSqScratch_->data()[I];
      DIFF = std::conj(DIFF)*DIFF;
      PARMS += reinterpret_cast<double(&)[2]>(DIFF)[0];
    }
  } else {
    for(auto I = 0; I < this->onePDMScalar_->size(); I++){
      T DIFF = this->onePDMScalar_->data()[I] - this->NBSqScratch_->data()[I];
      DIFF = std::conj(DIFF)*DIFF;
      PARMS += reinterpret_cast<double(&)[2]>(DIFF)[0];
    }

    DMzOld_->read(this->NBSqScratch_->data(),H5PredType<T>());
    for(auto I = 0; I < this->onePDMMz_->size(); I++){
      T DIFF = this->onePDMMz_->data()[I] - this->NBSqScratch_->data()[I];
      DIFF = std::conj(DIFF)*DIFF;
      PBRMS += reinterpret_cast<double(&)[2]>(DIFF)[0];
    }
  }


  PARMS = std::sqrt(PARMS);
  if(this->nTCS_ == 1 && !this->isClosedShell)
    PBRMS = std::sqrt(PBRMS);
  
  SCFConvergence CONVER;
  CONVER.EDelta = EDelta;
  CONVER.PARMS  = PARMS;
  if(this->nTCS_ == 1 && !this->isClosedShell)
    CONVER.PBRMS = PBRMS;
  
  EDelta = std::abs(EDelta);
  this->isConverged = EDelta < this->eneTol_;
  this->isConverged = this->isConverged && PARMS < this->denTol_;
  if(this->nTCS_ == 1 && !this->isClosedShell)
    this->isConverged = this->isConverged && PBRMS < this->denTol_;

  this->isConverged = this->isConverged || EDelta < this->eneTol_*1e-3;
  
  if(this->isPrimary) this->writeSCFFiles();

  return CONVER;
};
