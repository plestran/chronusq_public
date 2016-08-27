/**
 * \brief Initialize the pointers that only occur in the SCF to NULL
 */
template <typename T>
void SingleSlater<T>::initSCFPtr(){
  this->occNumMem_     = NULL;

  this->lambdaMem_     = NULL;
  this->delFMem_       = NULL;
  this->PNOMem_        = NULL;

}; // initSCFPtr

/**
 *  \brief Allocate RAM and Disc for the SCF Procedure
 */
template<typename T>
void SingleSlater<T>::initSCFMem3(){
  // Initialize the pointers to NULL
  this->initSCFPtr();

  auto NSQ = this->nBasis_ * this->nBasis_; 


  if(this->Ref_ == CUHF) {
    this->delFMem_   = this->memManager_->template malloc<T>(NSQ);
    this->lambdaMem_ = this->memManager_->template malloc<T>(NSQ);
    this->PNOMem_    = this->memManager_->template malloc<T>(NSQ);
    this->occNumMem_ = 
      this->memManager_->template malloc<double>(this->nBasis_);

    // DIIS NYI for CUHF
    this->doDIIS = false;
    cout << "IM HERE IN CUHF" << endl;
  }

  // Initalize the files for DIIS and DMS 
  if(this->doDIIS) this->initDIISFiles();
  if(this->doDMS)  this->initDMSFiles();

  // Allocate Files for SCF
  std::vector<hsize_t> dims;
  dims.push_back(this->nBasis_);
  dims.push_back(this->nBasis_);

  // Scalar Files
  this->DScalarOld_ = this->fileio_->createScratchPartition(H5PredType<T>(),
    "Most Recent Density Matrix (Scalar) for RMS",dims);
  this->PTScalarOld_ = this->fileio_->createScratchPartition(H5PredType<T>(),
    "Most Recent Perturbation Tensor (Scalar)",dims);
  this->DeltaDScalar_ = this->fileio_->createScratchPartition(H5PredType<T>(),
    "Change in Density Matrix (Scalar)",dims);

  // Mz Files
  if(this->nTCS_ == 2 || !this->isClosedShell){
    this->DMzOld_ = this->fileio_->createScratchPartition(H5PredType<T>(),
      "Most Recent Density Matrix (Mz) for RMS",dims);
    this->PTMzOld_ = this->fileio_->createScratchPartition(H5PredType<T>(),
      "Most Recent Perturbation Tensor (Mz)",dims);
    this->DeltaDMz_ = this->fileio_->createScratchPartition(H5PredType<T>(),
      "Change in Density Matrix (Mz)",dims);
  }
  
  // Mx My Files
  if(this->nTCS_ == 2) {
    this->DMyOld_ = this->fileio_->createScratchPartition(H5PredType<T>(),
      "Most Recent Density Matrix (My) for RMS",dims);
    this->DMxOld_ = this->fileio_->createScratchPartition(H5PredType<T>(),
      "Most Recent Density Matrix (Mx) for RMS",dims);
    this->PTMyOld_ = this->fileio_->createScratchPartition(H5PredType<T>(),
      "Most Recent Perturbation Tensor (My)",dims);
    this->PTMxOld_ = this->fileio_->createScratchPartition(H5PredType<T>(),
      "Most Recent Perturbation Tensor (Mx)",dims);
    this->DeltaDMy_ = this->fileio_->createScratchPartition(H5PredType<T>(),
      "Change in Density Matrix (My)",dims);
    this->DeltaDMx_ = this->fileio_->createScratchPartition(H5PredType<T>(),
      "Change in Density Matrix (Mx)",dims);
  }

};

/**
 * \brief Deallocate memory associated with SCF
 */
template<typename T>
void SingleSlater<T>::cleanupSCFMem3(){

  auto NSQ = this->nBasis_  * this->nBasis_;

  if(this->Ref_ == CUHF){ 
    this->memManager_->free(this->delFMem_,NSQ);   
    this->memManager_->free(this->lambdaMem_,NSQ); 
    this->memManager_->free(this->PNOMem_,NSQ);    
    this->memManager_->free(this->occNumMem_,this->nBasis_); 
  }
}

/**
 * \brief Copy the current quaternion density matrix to disc
 */
template <typename T>
void SingleSlater<T>::copyDen(){

  // Scalar
  DScalarOld_->write(this->onePDMScalar_->data(),H5PredType<T>());

  // Mz
  if(this->nTCS_ == 2 or !this->isClosedShell){
    DMzOld_->write(this->onePDMMz_->data(),H5PredType<T>());
  }

  // Mx My
  if(this->nTCS_ == 2){
    DMyOld_->write(this->onePDMMy_->data(),H5PredType<T>());
    DMxOld_->write(this->onePDMMx_->data(),H5PredType<T>());
  }
};

/**
 *  \brief Orthogonalize the quaternion Fock matrix 
 */
template<typename T>
void SingleSlater<T>::orthoFock3(){

  // Scalar
  this->aointegrals_->Ortho1Trans(*this->fockScalar_,*this->fockOrthoScalar_);
  
  if(this->nTCS_ == 2 or !this->isClosedShell) {
    // Mz
    this->aointegrals_->Ortho1Trans(*this->fockMz_,*this->fockOrthoMz_);
  }

  // Mx My
  if(this->nTCS_ == 2){
    this->aointegrals_->Ortho1Trans(*this->fockMy_,*this->fockOrthoMx_);
    this->aointegrals_->Ortho1Trans(*this->fockMy_,*this->fockOrthoMx_);
  }

  // Populate the gathered Fock matricies 
//this->gatherOrthoFock();

}

/**
 *  \brief Transform the quaternion orthonormal Fock matricies to
 *  their gathered form
 */
/*
template<typename T>
void SingleSlater<T>::gatherOrthoFock(){
  // FIXME: This logic needs to be fixed for RHF when RHF is just scalar 
  if(this->nTCS_ == 1 && this->isClosedShell) return;

  std::vector<std::reference_wrapper<TMap>> toGather;
  toGather.emplace_back(*this->fockOrthoScalar_);
  toGather.emplace_back(*this->fockOrthoMz_);
  if(this->nTCS_ == 1)
    Quantum<T>::spinGather(*this->fockOrthoA_,*this->fockOrthoB_,toGather);
  else {
    toGather.emplace_back(*this->fockOrthoMy_);
    toGather.emplace_back(*this->fockOrthoMx_);
    Quantum<T>::spinGather(*this->fockOrthoA_,toGather);
  }

};
*/

/**
 *  \brief Transform the orthonormal quaternion density to the AO basis
 */
template<typename T>
void SingleSlater<T>::unOrthoDen3(){
  this->aointegrals_->Ortho1TransT(
    *this->onePDMOrthoScalar_,*this->onePDMScalar_);

  if(this->nTCS_ == 2 or !this->isClosedShell){
    this->aointegrals_->Ortho1TransT(
      *this->onePDMOrthoMz_,*this->onePDMMz_);
  }

  if(this->nTCS_ == 2){
    this->aointegrals_->Ortho1TransT(
      *this->onePDMOrthoMy_,*this->onePDMMx_);
    this->aointegrals_->Ortho1TransT(
      *this->onePDMOrthoMy_,*this->onePDMMx_);
  }
}


/**
 *  \brief Copy the AO quaternion density matrix to the orthonormal
 *  quaternion storage
 *
 *  This is a hack for the SCF as we build the density in the ortho
 *  normal basis to avoid transforming the full 2C MO coeffients at
 *  every SCF iteration. Because formDensity populates the AO 2C and
 *  quaternion storage, it must be copied into it's correct place 
 *  after the call in SCF
 */
template<typename T>
void SingleSlater<T>::cpyAOtoOrthoDen(){
/*
  if(this->nTCS_ == 1 && this->isClosedShell) {
    // Just copy for RHF
    (*this->onePDMOrthoScalar_) = (*this->onePDMScalar_);
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
*/
  (*this->onePDMOrthoScalar_) = (*this->onePDMScalar_);
  if(this->nTCS_ == 2 or !this->isClosedShell)
    (*this->onePDMOrthoMz_) = (*this->onePDMMz_);
  if(this->nTCS_ == 2) {
    (*this->onePDMOrthoMy_) = (*this->onePDMMy_);
    (*this->onePDMOrthoMx_) = (*this->onePDMMx_);
  }

};

/**
 * \brief Evaluate the convergence criteria for the SCF procedure
 */
template<typename T>
SCFConvergence SingleSlater<T>::evalConver3(){

  // Energy Convergence
  double EOld = this->totalEnergy;
  this->computeEnergy();
  double EDelta = this->totalEnergy - EOld;

  double PSRMS(0),PMRMS(0);
  // Write D(M) - D(M-1) to disc
  this->formDeltaD();

  // Scalar density RMS difference
  DeltaDScalar_->read(this->NBSqScratch_->data(),H5PredType<T>());
  PSRMS = this->NBSqScratch_->norm();

  // Magnetization density RMS
     
  // || Pz(M) - Pz(M-1) ||^2
  if(this->nTCS_ == 2 or !this->isClosedShell){
    DeltaDMz_->read(this->NBSqScratch_->data(),H5PredType<T>());
    PMRMS = this->NBSqScratch_->squaredNorm();
  }

  // || Py(M) - Py(M-1) ||^2 + || Px(M) - Px(M-1) ||^2
  if(this->nTCS_ == 2) {
    DeltaDMy_->read(this->NBSqScratch_->data(),H5PredType<T>());
    PMRMS += this->NBSqScratch_->squaredNorm();
    DeltaDMx_->read(this->NBSqScratch_->data(),H5PredType<T>());
    PMRMS += this->NBSqScratch_->squaredNorm();
  }

  // Sqrt of norm squared
  PMRMS = std::sqrt(PMRMS);

  
  // Check if Density and Energy are properly converged
  SCFConvergence CONVER;
  CONVER.EDelta = EDelta;
  CONVER.PSRMS  = PSRMS;
  if(this->nTCS_ == 1 && !this->isClosedShell)
    CONVER.PMRMS = PMRMS;
  
  if(!this->isConverged) { // [F,P] supercedes this check
    EDelta = std::abs(EDelta);
    this->isConverged = EDelta < this->eneTol_;
    this->isConverged = this->isConverged && PSRMS < this->denTol_;
    if(this->nTCS_ == 2 or !this->isClosedShell)
      this->isConverged = this->isConverged && PMRMS < this->denTol_;
  }

//this->isConverged = this->isConverged || EDelta < this->eneTol_*1e-3;
  
  // Copy the current SingleSlater moieties to disc (this is wrong for now)
//if(this->isPrimary) this->writeSCFFiles();

  return CONVER;
};

/**
 * \brief Properly mix orbitals for Complex / 2C guess
 */
template<typename T>
void SingleSlater<T>::mixOrbitalsSCF(){
  this->mixOrbitals2C();
  this->mixOrbitalsComplex();
}

/**
 *  \brief Properly mix orbitals for 2C guess
 */
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

/**
 * \brief Write the change in the quaternion density to Disc 
 */
template <typename T>
void SingleSlater<T>::formDeltaD(){
  DScalarOld_->read(this->NBSqScratch_->data(),H5PredType<T>());
  (*this->NBSqScratch2_) = (*this->onePDMScalar_) - (*this->NBSqScratch_);
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

/**
 * \brief Copy the Most recent change in density to the quaternion density
 * storage
 *
 * This is for use with the incremental Fock build to populate the quaternion
 * density to be traced with the ERIs to compute the change of F at the current
 * SCF iteration
 */
template<typename T>
void SingleSlater<T>::copyDeltaDtoD(){
  DeltaDScalar_->read(this->onePDMScalar_->data(),H5PredType<T>());
  if(this->nTCS_ == 2 or !this->isClosedShell){
    DeltaDMz_->read(this->onePDMMz_->data(),H5PredType<T>());
  }
  if(this->nTCS_ == 2) {
    DeltaDMy_->read(this->onePDMMy_->data(),H5PredType<T>());
    DeltaDMx_->read(this->onePDMMx_->data(),H5PredType<T>());
  }
}

/**
 * \brief Repopulate the quaternion density matrix with the most recent
 * density written to disk
 *
 * This is for use with the incremental Fock build to repopulate the quaternion
 * density after the change in F has been evaluated 
 */
template<typename T>
void SingleSlater<T>::copyDOldtoD(){
  DScalarOld_->read(this->onePDMScalar_->data(),H5PredType<T>());
  if(this->nTCS_ == 2 or !this->isClosedShell){
    DMzOld_->read(this->onePDMMz_->data(),H5PredType<T>());
  }
  if(this->nTCS_ == 2) {
    DMyOld_->read(this->onePDMMy_->data(),H5PredType<T>());
    DMxOld_->read(this->onePDMMx_->data(),H5PredType<T>());
  }
}

/**
 * \brief Write the most recent quaternion perturbation tensor to disk
 *
 * This is for use with the incremental Fock build to compute the total
 * quaternion perturbation tensor after only the change in the perturbation
 * tensor has been evalued in formFock
 */
template <typename T>
void SingleSlater<T>::copyPT(){

  PTScalarOld_->write(this->PTScalar_->data(),H5PredType<T>());
  if(this->nTCS_ == 2 or !this->isClosedShell){
    PTMzOld_->write(this->PTMz_->data(),H5PredType<T>());
  }
  if(this->nTCS_ == 2){
    PTMyOld_->write(this->PTMy_->data(),H5PredType<T>());
    PTMxOld_->write(this->PTMx_->data(),H5PredType<T>());
  }
};

/**
 *  \brief Increment the current quaternion perturbation tensor by the last
 *  quaternion perturbation tensor written to disk
 *
 *  This is for use with the incremental Fock build to repopulate the full
 *  quaternion perturbation tensor
 */
template<typename T>
void SingleSlater<T>::incPT(){
  PTScalarOld_->read(this->NBSqScratch_->data(),H5PredType<T>());
  this->PTScalar_->noalias() += (*this->NBSqScratch_);

  if(this->nTCS_ == 2 or !this->isClosedShell){
    PTMzOld_->read(this->NBSqScratch_->data(),H5PredType<T>());
    this->PTMz_->noalias() += (*this->NBSqScratch_);
  }

  if(this->nTCS_ == 2) {
    PTMxOld_->read(this->NBSqScratch_->data(),H5PredType<T>());
    this->PTMx_->noalias() += (*this->NBSqScratch_);
    PTMyOld_->read(this->NBSqScratch_->data(),H5PredType<T>());
    this->PTMy_->noalias() += (*this->NBSqScratch_);
  }
}

template<typename T>
void SingleSlater<T>::populateMO4Diag(){
  if(this->nTCS_ == 1 and this->isClosedShell)
    (*this->moA_) = 0.5 * (*this->fockOrthoScalar_);
  else if(this->nTCS_ == 1 and !this->isClosedShell) {
    (*this->moA_) = 0.5 * ((*this->fockOrthoScalar_) + (*this->fockOrthoMz_));
    (*this->moB_) = 0.5 * ((*this->fockOrthoScalar_) - (*this->fockOrthoMz_));
  } else {
    std::vector<std::reference_wrapper<TMap>> scattered;
    for(auto iF = this->fockOrtho_.begin(); iF != this->fockOrtho_.end(); iF++)
      scattered.emplace_back(*(*iF));
    Quantum<T>::spinGather(*this->moA_,scattered);
  }
};
