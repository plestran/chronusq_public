template <typename T>
void SingleSlater<T>::initSCFPtr(){
  this->occNumMem_     = NULL;

/*
  this->FpAlphaMem_    = NULL;
  this->FpBetaMem_     = NULL;
  this->POldAlphaMem_  = NULL;
  this->POldBetaMem_   = NULL;
  this->ErrorAlphaMem_ = NULL;
  this->ErrorBetaMem_  = NULL;
  this->FADIIS_        = NULL;
  this->FBDIIS_        = NULL;
*/
  this->lambdaMem_     = NULL;
  this->delFMem_       = NULL;
  this->PNOMem_        = NULL;

/*
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
*/
}; // initSCFPtr

template<typename T>
void SingleSlater<T>::initSCFMem3(){
  this->initSCFPtr();

  auto NTCSxNBASIS = this->nTCS_ * this->nBasis_;
  auto NSQ = NTCSxNBASIS  * NTCSxNBASIS;

/*
  this->POldAlphaMem_  = this->memManager_->template malloc<T>(NSQ);
  if(this->nTCS_ == 1 && !this->isClosedShell) 
    this->POldBetaMem_  = this->memManager_->template malloc<T>(NSQ);
*/

  if(this->Ref_ == CUHF) {
    this->delFMem_   = this->memManager_->template malloc<T>(NSQ);
    this->lambdaMem_ = this->memManager_->template malloc<T>(NSQ);
    this->PNOMem_    = this->memManager_->template malloc<T>(NSQ);
    this->occNumMem_ = this->memManager_->template malloc<double>(NTCSxNBASIS);
    this->doDIIS = false;
    cout << "IM HERE IN CUHF" << endl;
  }

   
  if(this->doDIIS) this->initDIISFiles();
  if(this->doDMS)  this->initDMSFiles();

  std::vector<hsize_t> dims;
  dims.push_back(this->nBasis_);
  dims.push_back(this->nBasis_);

  this->DScalarOld_ = this->fileio_->createScratchPartition(H5PredType<T>(),
    "Most Recent Density Matrix (Scalar) for RMS",dims);
  this->DeltaDScalar_ = this->fileio_->createScratchPartition(H5PredType<T>(),
    "Change in Density Matrix (Scalar)",dims);

  if(this->nTCS_ == 2 || !this->isClosedShell){
    this->DMzOld_ = this->fileio_->createScratchPartition(H5PredType<T>(),
      "Most Recent Density Matrix (Mz) for RMS",dims);
    this->DeltaDMz_ = this->fileio_->createScratchPartition(H5PredType<T>(),
      "Change in Density Matrix (Mz)",dims);
  }
  
  if(this->nTCS_ == 2) {
    this->DMyOld_ = this->fileio_->createScratchPartition(H5PredType<T>(),
      "Most Recent Density Matrix (My) for RMS",dims);
    this->DMxOld_ = this->fileio_->createScratchPartition(H5PredType<T>(),
      "Most Recent Density Matrix (Mx) for RMS",dims);
    this->DeltaDMy_ = this->fileio_->createScratchPartition(H5PredType<T>(),
      "Change in Density Matrix (My)",dims);
    this->DeltaDMx_ = this->fileio_->createScratchPartition(H5PredType<T>(),
      "Change in Density Matrix (Mx)",dims);
  }

};

template<typename T>
void SingleSlater<T>::cleanupSCFMem3(){

  auto NTCSxNBASIS = this->nTCS_ * this->nBasis_;
  auto NSQ = NTCSxNBASIS  * NTCSxNBASIS;
/*
  this->memManager_->free(this->POldAlphaMem_,NSQ);  
  if(this->nTCS_ == 1 && !this->isClosedShell) 
    this->memManager_->free(this->POldBetaMem_,NSQ);  
*/

  if(this->Ref_ == CUHF){ 
    this->memManager_->free(this->delFMem_,NSQ);   
    this->memManager_->free(this->lambdaMem_,NSQ); 
    this->memManager_->free(this->PNOMem_,NSQ);    
    this->memManager_->free(this->occNumMem_,NTCSxNBASIS); 
  }
}

template <typename T>
void SingleSlater<T>::copyDen(){
  T* ScalarPtr;
  if(this->isClosedShell && this->nTCS_ == 1) 
    ScalarPtr = this->onePDMA_->data();
  else
    ScalarPtr = this->onePDMScalar_->data();

  DScalarOld_->write(ScalarPtr,H5PredType<T>());
  if(this->nTCS_ == 2 || !this->isClosedShell){
    DMzOld_->write(this->onePDMMz_->data(),H5PredType<T>());
  }
  if(this->nTCS_ == 2){
    DMyOld_->write(this->onePDMMy_->data(),H5PredType<T>());
    DMxOld_->write(this->onePDMMx_->data(),H5PredType<T>());
  }
};

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
    this->aointegrals_->Ortho1TransT(
      *this->onePDMOrthoA_,*this->onePDMA_);
  } else {
    this->aointegrals_->Ortho1TransT(
      *this->onePDMOrthoScalar_,*this->onePDMScalar_);
    this->aointegrals_->Ortho1TransT(
      *this->onePDMOrthoMz_,*this->onePDMMz_);
    if(this->nTCS_ == 2){
      this->aointegrals_->Ortho1TransT(
        *this->onePDMOrthoMy_,*this->onePDMMx_);
      this->aointegrals_->Ortho1TransT(
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
  }
};

template<typename T>
SCFConvergence SingleSlater<T>::evalConver3(){

  // Energy Convergence
  double EOld = this->totalEnergy;
//this->formFock();
  this->computeEnergy();
  double EDelta = this->totalEnergy - EOld;

  double PARMS(0),PBRMS(0);
  this->formDeltaD();

/*
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
*/

  
  DeltaDScalar_->read(this->NBSqScratch_->data(),H5PredType<T>());
  PARMS = this->NBSqScratch_->norm();
  if(this->nTCS_ == 2 or !this->isClosedShell){
    DeltaDMz_->read(this->NBSqScratch_->data(),H5PredType<T>());
    PBRMS = this->NBSqScratch_->squaredNorm();
  }
  if(this->nTCS_ == 2) {
    DeltaDMy_->read(this->NBSqScratch_->data(),H5PredType<T>());
    PBRMS += this->NBSqScratch_->squaredNorm();
    DeltaDMx_->read(this->NBSqScratch_->data(),H5PredType<T>());
    PBRMS += this->NBSqScratch_->squaredNorm();
  }
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

//this->isConverged = this->isConverged || EDelta < this->eneTol_*1e-3;
  
  if(this->isPrimary) this->writeSCFFiles();

  return CONVER;
};

template<typename T>
void SingleSlater<T>::mixOrbitalsSCF(){
  this->mixOrbitals2C();
  this->mixOrbitalsComplex();
}

template<typename T>
void SingleSlater<T>::mixOrbitals2C(){
  auto nO = this->nAE_ + this->nBE_;
  if(this->nTCS_ != 2) return;
  this->fileio_->out << "** Mixing Alpha-Beta Orbitals for 2C Guess **" << endl;
  TVec HOMOA,LUMOB;
  int indxHOMOA = -1, indxLUMOB = -1;
  auto nOrb = this->nBasis_;
  double maxPercentNonZeroAlpha = 0;

  for(auto i = nO-1; i >= 0; i--){
    auto nNonZeroAlpha = 0;
    for(auto j = 0; j < this->nTCS_*this->nBasis_; j+=2){
      auto aComp = (*this->moA_)(j,i);
      auto bComp = (*this->moA_)(j+1,i);
      if(std::norm(aComp) > 1e-12 && std::norm(bComp) < 1e-12) nNonZeroAlpha++;
    }
    double percentNonZeroAlpha = (double)nNonZeroAlpha/(double)nOrb;
    if(percentNonZeroAlpha > maxPercentNonZeroAlpha){
      maxPercentNonZeroAlpha = percentNonZeroAlpha;
      indxHOMOA = i;
    }
  }

  double maxPercentNonZeroBeta = 0;
  for(auto i = nO; i < this->nTCS_*this->nBasis_; i++){
    auto nNonZeroBeta = 0;
    for(auto j = 1; j < this->nTCS_*this->nBasis_; j+=2){
      auto aComp = (*this->moA_)(j-1,i);
      auto bComp = (*this->moA_)(j,i);
      if(std::norm(bComp) > 1e-12 && std::norm(aComp) < 1e-12) nNonZeroBeta++;
    }
    double percentNonZeroBeta = (double)nNonZeroBeta/(double)nOrb;
    if(percentNonZeroBeta > maxPercentNonZeroBeta){
      maxPercentNonZeroBeta = percentNonZeroBeta;
      indxLUMOB = i;
    }
  }

  if(indxHOMOA == -1 || indxLUMOB == -1){
    this->fileio_->out 
      << "TCS orbital swap failed to find suitable Alpha-Beta pair" << endl;
    return;
  }
  
  HOMOA = this->moA_->col(indxHOMOA) ;
  LUMOB = this->moA_->col(indxLUMOB) ;
  this->moA_->col(indxHOMOA) = std::sqrt(0.5) * (HOMOA + LUMOB);
  this->moA_->col(indxLUMOB) = std::sqrt(0.5) * (HOMOA - LUMOB);
}

template <typename T>
void SingleSlater<T>::formDeltaD(){
  DScalarOld_->read(this->NBSqScratch_->data(),H5PredType<T>());
  if(this->nTCS_ == 1 && this->isClosedShell){
    (*this->NBSqScratch2_) = (*this->onePDMA_) - (*this->NBSqScratch_);
  } else {
    (*this->NBSqScratch2_) = (*this->onePDMScalar_) - (*this->NBSqScratch_);
  }

  DeltaDScalar_->write(this->NBSqScratch2_->data(),H5PredType<T>());

  if(this->nTCS_ == 2 or !this->isClosedShell){
    DMzOld_->read(this->NBSqScratch_->data(),H5PredType<T>());
    (*this->NBSqScratch2_) = (*this->onePDMMz_) - (*this->NBSqScratch_);
    DeltaDMz_->write(this->NBSqScratch2_->data(),H5PredType<T>());
  }

  if(this->nTCS_ == 2) {
    DMyOld_->read(this->NBSqScratch_->data(),H5PredType<T>());
    (*this->NBSqScratch2_) = (*this->onePDMMy_) - (*this->NBSqScratch_);
    DeltaDMy_->write(this->NBSqScratch2_->data(),H5PredType<T>());

    DMxOld_->read(this->NBSqScratch_->data(),H5PredType<T>());
    (*this->NBSqScratch2_) = (*this->onePDMMx_) - (*this->NBSqScratch_);
    DeltaDMx_->write(this->NBSqScratch2_->data(),H5PredType<T>());
  }

};
