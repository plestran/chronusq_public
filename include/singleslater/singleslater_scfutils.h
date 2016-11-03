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
//this->SCFDensityScalar_ = this->fileio_->createScratchPartition(H5PredType<T>(),
//  "Most Recent Density Matrix (Scalar) for RMS",dims);
//this->PTScalarOld_ = this->fileio_->createScratchPartition(H5PredType<T>(),
//  "Most Recent Perturbation Tensor (Scalar)",dims);
  this->DeltaDScalar_ = this->fileio_->createScratchPartition(H5PredType<T>(),
    "Change in Density Matrix (Scalar)",dims);

  // Mz Files
  if(this->nTCS_ == 2 || !this->isClosedShell){
//  this->SCFDensityMz_ = this->fileio_->createScratchPartition(H5PredType<T>(),
//    "Most Recent Density Matrix (Mz) for RMS",dims);
//  this->PTMzOld_ = this->fileio_->createScratchPartition(H5PredType<T>(),
//    "Most Recent Perturbation Tensor (Mz)",dims);
    this->DeltaDMz_ = this->fileio_->createScratchPartition(H5PredType<T>(),
      "Change in Density Matrix (Mz)",dims);
  }
  
  // Mx My Files
  if(this->nTCS_ == 2) {
//  this->SCFDensityMy_ = this->fileio_->createScratchPartition(H5PredType<T>(),
//    "Most Recent Density Matrix (My) for RMS",dims);
//  this->SCFDensityMx_ = this->fileio_->createScratchPartition(H5PredType<T>(),
//    "Most Recent Density Matrix (Mx) for RMS",dims);
//  this->PTMyOld_ = this->fileio_->createScratchPartition(H5PredType<T>(),
//    "Most Recent Perturbation Tensor (My)",dims);
//  this->PTMxOld_ = this->fileio_->createScratchPartition(H5PredType<T>(),
//    "Most Recent Perturbation Tensor (Mx)",dims);
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
/*
template <typename T>
void SingleSlater<T>::copyDen(){

  // Scalar
  SCFDensityScalar_->write(this->onePDMScalar_->data(),H5PredType<T>());

  // Mz
  if(this->nTCS_ == 2 or !this->isClosedShell){
    SCFDensityMz_->write(this->onePDMMz_->data(),H5PredType<T>());
  }

  // Mx My
  if(this->nTCS_ == 2){
    SCFDensityMy_->write(this->onePDMMy_->data(),H5PredType<T>());
    SCFDensityMx_->write(this->onePDMMx_->data(),H5PredType<T>());
  }
};
*/

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
    this->aointegrals_->Ortho1Trans(*this->fockMy_,*this->fockOrthoMy_);
    this->aointegrals_->Ortho1Trans(*this->fockMx_,*this->fockOrthoMx_);
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
      *this->onePDMOrthoMx_,*this->onePDMMx_);
    this->aointegrals_->Ortho1TransT(
      *this->onePDMOrthoMy_,*this->onePDMMy_);
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
  double EOld = this->totalEnergy_;
  this->computeEnergy();
  double EDelta = this->totalEnergy_ - EOld;

  // Turn off damping when close to convergence
  if(this->doDamp and (std::abs(EDelta) < 1e-6)){
    if(this->dampParam != 0.)
      this->fileio_->out << 
        "    *** Damping Disabled After 1e-6 Converged Met ***" << endl;
    this->dampParam = 0.0;
  } else if(this->doDamp and (std::abs(EDelta) > 1e-6) and 
            this->dampParam <= 0.){
    this->fileio_->out << 
      "    *** Damping Enabled due to > 1e-6 Oscillation in Energy ***" << endl;
    this->dampParam = 0.1;
  }

  // Turn off level shifting when close to convergence
  if(this->doLevelShift and (std::abs(EDelta) < 1e-5)){
    if(this->levelShiftParam != 0.)
      this->fileio_->out << 
        "    *** VShift Disabled After 1e-6 Converged Met ***" << endl;
    this->levelShiftParam = 0.0;
  } else if(this->doLevelShift and (std::abs(EDelta) > 1e-5) and 
            this->levelShiftParam <= 0.){
    this->fileio_->out << 
      "    *** VShift Enabled due to > 1e-6 Oscillation in Energy ***" << endl;
    this->levelShiftParam = 0.1;
  }

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
  if(this->nTCS_ != 2) return;
  this->fileio_->out << "** Mixing Alpha-Beta Orbitals for 2C Guess **" << endl;
  TVec HOMOA,LUMOB;
  int indxHOMOA = -1, indxLUMOB = -1;
  auto nOrb = this->nBasis_;
  double maxPercentNonZeroAlpha = 0;

  for(auto i = this->nO_ - 1; i >= 0; i--){
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
  for(auto i = this->nO_; i < this->nTCS_*this->nBasis_; i++){
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
  SCFDensityScalar_->read(this->NBSqScratch_->data(),H5PredType<T>());
  (*this->NBSqScratch2_) = (*this->onePDMScalar_) - (*this->NBSqScratch_);
  DeltaDScalar_->write(this->NBSqScratch2_->data(),H5PredType<T>());

  if(this->nTCS_ == 2 or !this->isClosedShell){
    SCFDensityMz_->read(this->NBSqScratch_->data(),H5PredType<T>());
    (*this->NBSqScratch2_) = (*this->onePDMMz_) - (*this->NBSqScratch_);
    DeltaDMz_->write(this->NBSqScratch2_->data(),H5PredType<T>());
  }

  if(this->nTCS_ == 2) {
    SCFDensityMy_->read(this->NBSqScratch_->data(),H5PredType<T>());
    (*this->NBSqScratch2_) = (*this->onePDMMy_) - (*this->NBSqScratch_);
    DeltaDMy_->write(this->NBSqScratch2_->data(),H5PredType<T>());

    SCFDensityMx_->read(this->NBSqScratch_->data(),H5PredType<T>());
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
  SCFDensityScalar_->read(this->onePDMScalar_->data(),H5PredType<T>());
  if(this->nTCS_ == 2 or !this->isClosedShell){
    SCFDensityMz_->read(this->onePDMMz_->data(),H5PredType<T>());
  }
  if(this->nTCS_ == 2) {
    SCFDensityMy_->read(this->onePDMMy_->data(),H5PredType<T>());
    SCFDensityMx_->read(this->onePDMMx_->data(),H5PredType<T>());
  }
}

/**
 * \brief Write the most recent quaternion perturbation tensor to disk
 *
 * This is for use with the incremental Fock build to compute the total
 * quaternion perturbation tensor after only the change in the perturbation
 * tensor has been evalued in formFock
 */
/*
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
*/

/**
 *  \brief Increment the current quaternion perturbation tensor by the last
 *  quaternion perturbation tensor written to disk
 *
 *  This is for use with the incremental Fock build to repopulate the full
 *  quaternion perturbation tensor
 */
template<typename T>
void SingleSlater<T>::incPT(){
  SCFPTScalar_->read(this->NBSqScratch_->data(),H5PredType<T>());
  this->PTScalar_->noalias() += (*this->NBSqScratch_);

  if(this->nTCS_ == 2 or !this->isClosedShell){
    SCFPTMz_->read(this->NBSqScratch_->data(),H5PredType<T>());
    this->PTMz_->noalias() += (*this->NBSqScratch_);
  }

  if(this->nTCS_ == 2) {
    SCFPTMx_->read(this->NBSqScratch_->data(),H5PredType<T>());
    this->PTMx_->noalias() += (*this->NBSqScratch_);
    SCFPTMy_->read(this->NBSqScratch_->data(),H5PredType<T>());
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
    int I = 0;
    for(auto iF = this->fockOrtho_.begin(); iF != this->fockOrtho_.end(); iF++){
      scattered.emplace_back(*(*iF));
    }
    Quantum<T>::spinGather(*this->moA_,scattered);
  }
};

template <typename T>
void SingleSlater<T>::initSCFFiles() {

  if(!this->isPrimary) return;
  std::vector<hsize_t> dims;
  dims.push_back(this->basisset_->nBasis());
  dims.push_back(this->basisset_->nBasis());

  std::vector<hsize_t> bigDims;
  bigDims.push_back(this->nTCS_*this->basisset_->nBasis());
  bigDims.push_back(this->nTCS_*this->basisset_->nBasis());

  H5::DataSpace dsp(2,&dims[0]);
  H5::DataSpace bdsp(2,&bigDims[0]);

  // Check if SCF group exists
  try {
    this->SCFGroup_ = std::unique_ptr<H5::Group>(new H5::Group(
      this->fileio_->restart->openGroup("/SCF")));
  } catch(...) {
    this->SCFGroup_ = std::unique_ptr<H5::Group>(new H5::Group(
      this->fileio_->restart->createGroup("/SCF")));
  }


  // Density Files

  // Scalar
  try {
    this->SCFDensityScalar_ = std::unique_ptr<H5::DataSet>( new H5::DataSet(
      this->SCFGroup_->openDataSet("DensityScalar")));
  } catch(...) {
    if(this->guess_ == READ) 
      CErr("FATAL: READ Guess cannot find DensityScalar",this->fileio_->out);

    this->SCFDensityScalar_ = std::unique_ptr<H5::DataSet>( new H5::DataSet(
      this->SCFGroup_->createDataSet("DensityScalar",H5PredType<T>(),dsp)));
  }

  // Mz
  if(this->nTCS_ == 2 or !this->isClosedShell) {
    try {
      this->SCFDensityMz_ = std::unique_ptr<H5::DataSet>( new H5::DataSet(
        this->SCFGroup_->openDataSet("DensityMz")));
    } catch(...) {
      this->SCFDensityMz_ = std::unique_ptr<H5::DataSet>( new H5::DataSet(
        this->SCFGroup_->createDataSet("DensityMz",H5PredType<T>(),dsp)));
    }
  }

  // Mx / My
  if(this->nTCS_ == 2) {
    // Mx
    try {
      this->SCFDensityMx_ = std::unique_ptr<H5::DataSet>( new H5::DataSet(
        this->SCFGroup_->openDataSet("DensityMx")));
    } catch(...) {
      this->SCFDensityMx_ = std::unique_ptr<H5::DataSet>( new H5::DataSet(
        this->SCFGroup_->createDataSet("DensityMx",H5PredType<T>(),dsp)));
    }

    // My
    try {
      this->SCFDensityMy_ = std::unique_ptr<H5::DataSet>( new H5::DataSet(
        this->SCFGroup_->openDataSet("DensityMy")));
    } catch(...) {
      this->SCFDensityMy_ = std::unique_ptr<H5::DataSet>( new H5::DataSet(
        this->SCFGroup_->createDataSet("DensityMy",H5PredType<T>(),dsp)));
    }
  }

  // Ortho Density Files

  // Scalar
  try {
    this->SCFOrthoDScalar_ = std::unique_ptr<H5::DataSet>( new H5::DataSet(
      this->SCFGroup_->openDataSet("OrthoDScalar")));
  } catch(...) {
    this->SCFOrthoDScalar_ = std::unique_ptr<H5::DataSet>( new H5::DataSet(
      this->SCFGroup_->createDataSet("OrthoDScalar",H5PredType<T>(),dsp)));
  }

  // Mz
  if(this->nTCS_ == 2 or !this->isClosedShell) {
    try {
      this->SCFOrthoDMz_ = std::unique_ptr<H5::DataSet>( new H5::DataSet(
        this->SCFGroup_->openDataSet("OrthoDMz")));
    } catch(...) {
      this->SCFOrthoDMz_ = std::unique_ptr<H5::DataSet>( new H5::DataSet(
        this->SCFGroup_->createDataSet("OrthoDMz",H5PredType<T>(),dsp)));
    }
  }

  // Mx / My
  if(this->nTCS_ == 2) {
    // Mx
    try {
      this->SCFOrthoDMx_ = std::unique_ptr<H5::DataSet>( new H5::DataSet(
        this->SCFGroup_->openDataSet("OrthoDMx")));
    } catch(...) {
      this->SCFOrthoDMx_ = std::unique_ptr<H5::DataSet>( new H5::DataSet(
        this->SCFGroup_->createDataSet("OrthoDMx",H5PredType<T>(),dsp)));
    }

    // My
    try {
      this->SCFOrthoDMy_ = std::unique_ptr<H5::DataSet>( new H5::DataSet(
        this->SCFGroup_->openDataSet("OrthoDMy")));
    } catch(...) {
      this->SCFOrthoDMy_ = std::unique_ptr<H5::DataSet>( new H5::DataSet(
        this->SCFGroup_->createDataSet("OrthoDMy",H5PredType<T>(),dsp)));
    }
  }



  // Fock Files

  // Scalar
  try {
    this->SCFFockScalar_ = std::unique_ptr<H5::DataSet>( new H5::DataSet(
      this->SCFGroup_->openDataSet("FockScalar")));
  } catch(...) {
    this->SCFFockScalar_ = std::unique_ptr<H5::DataSet>( new H5::DataSet(
      this->SCFGroup_->createDataSet("FockScalar",H5PredType<T>(),dsp)));
  }

  // Mz
  if(this->nTCS_ == 2 or !this->isClosedShell) {
    try {
      this->SCFFockMz_ = std::unique_ptr<H5::DataSet>( new H5::DataSet(
        this->SCFGroup_->openDataSet("FockMz")));
    } catch(...) {
      this->SCFFockMz_ = std::unique_ptr<H5::DataSet>( new H5::DataSet(
        this->SCFGroup_->createDataSet("FockMz",H5PredType<T>(),dsp)));
    }
  }

  // Mx / My
  if(this->nTCS_ == 2) {
    // Mx
    try {
      this->SCFFockMx_ = std::unique_ptr<H5::DataSet>( new H5::DataSet(
        this->SCFGroup_->openDataSet("FockMx")));
    } catch(...) {
      this->SCFFockMx_ = std::unique_ptr<H5::DataSet>( new H5::DataSet(
        this->SCFGroup_->createDataSet("FockMx",H5PredType<T>(),dsp)));
    }

    // My
    try {
      this->SCFFockMy_ = std::unique_ptr<H5::DataSet>( new H5::DataSet(
        this->SCFGroup_->openDataSet("FockMy")));
    } catch(...) {
      this->SCFFockMy_ = std::unique_ptr<H5::DataSet>( new H5::DataSet(
        this->SCFGroup_->createDataSet("FockMy",H5PredType<T>(),dsp)));
    }
  }

  // Ortho Fock Files

  // Scalar
  try {
    this->SCFOrthoFScalar_ = std::unique_ptr<H5::DataSet>( new H5::DataSet(
      this->SCFGroup_->openDataSet("OrthoFScalar")));
  } catch(...) {
    this->SCFOrthoFScalar_ = std::unique_ptr<H5::DataSet>( new H5::DataSet(
      this->SCFGroup_->createDataSet("OrthoFScalar",H5PredType<T>(),dsp)));
  }

  // Mz
  if(this->nTCS_ == 2 or !this->isClosedShell) {
    try {
      this->SCFOrthoFMz_ = std::unique_ptr<H5::DataSet>( new H5::DataSet(
        this->SCFGroup_->openDataSet("OrthoFMz")));
    } catch(...) {
      this->SCFOrthoFMz_ = std::unique_ptr<H5::DataSet>( new H5::DataSet(
        this->SCFGroup_->createDataSet("OrthoFMz",H5PredType<T>(),dsp)));
    }
  }

  // Mx / My
  if(this->nTCS_ == 2) {
    // Mx
    try {
      this->SCFOrthoFMx_ = std::unique_ptr<H5::DataSet>( new H5::DataSet(
        this->SCFGroup_->openDataSet("OrthoFMx")));
    } catch(...) {
      this->SCFOrthoFMx_ = std::unique_ptr<H5::DataSet>( new H5::DataSet(
        this->SCFGroup_->createDataSet("OrthoFMx",H5PredType<T>(),dsp)));
    }

    // My
    try {
      this->SCFOrthoFMy_ = std::unique_ptr<H5::DataSet>( new H5::DataSet(
        this->SCFGroup_->openDataSet("OrthoFMy")));
    } catch(...) {
      this->SCFOrthoFMy_ = std::unique_ptr<H5::DataSet>( new H5::DataSet(
        this->SCFGroup_->createDataSet("OrthoFMy",H5PredType<T>(),dsp)));
    }
  }

  // PT Files

  // Scalar
  try {
    this->SCFPTScalar_ = std::unique_ptr<H5::DataSet>( new H5::DataSet(
      this->SCFGroup_->openDataSet("PTScalar")));
  } catch(...) {
    this->SCFPTScalar_ = std::unique_ptr<H5::DataSet>( new H5::DataSet(
      this->SCFGroup_->createDataSet("PTScalar",H5PredType<T>(),dsp)));
  }

  // Mz
  if(this->nTCS_ == 2 or !this->isClosedShell) {
    try {
      this->SCFPTMz_ = std::unique_ptr<H5::DataSet>( new H5::DataSet(
        this->SCFGroup_->openDataSet("PTMz")));
    } catch(...) {
      this->SCFPTMz_ = std::unique_ptr<H5::DataSet>( new H5::DataSet(
        this->SCFGroup_->createDataSet("PTMz",H5PredType<T>(),dsp)));
    }
  }

  // Mx / My
  if(this->nTCS_ == 2) {
    // Mx
    try {
      this->SCFPTMx_ = std::unique_ptr<H5::DataSet>( new H5::DataSet(
        this->SCFGroup_->openDataSet("PTMx")));
    } catch(...) {
      this->SCFPTMx_ = std::unique_ptr<H5::DataSet>( new H5::DataSet(
        this->SCFGroup_->createDataSet("PTMx",H5PredType<T>(),dsp)));
    }

    // My
    try {
      this->SCFPTMy_ = std::unique_ptr<H5::DataSet>( new H5::DataSet(
        this->SCFGroup_->openDataSet("PTMy")));
    } catch(...) {
      this->SCFPTMy_ = std::unique_ptr<H5::DataSet>( new H5::DataSet(
        this->SCFGroup_->createDataSet("PTMy",H5PredType<T>(),dsp)));
    }
  }

  // MO Files
  
  // MO (TCS) / MOA
  try {
    this->SCFMOA_ = std::unique_ptr<H5::DataSet>(new H5::DataSet(
      this->SCFGroup_->openDataSet("MOA")));
  } catch(...) {
    this->SCFMOA_ = std::unique_ptr<H5::DataSet>(new H5::DataSet(
      this->SCFGroup_->createDataSet("MOA",H5PredType<T>(),bdsp)));
  }

  // MOB
  if(this->nTCS_ == 1 and !this->isClosedShell){
    try {
      this->SCFMOB_ = std::unique_ptr<H5::DataSet>(new H5::DataSet(
        this->SCFGroup_->openDataSet("MOB")));
    } catch(...) {
      this->SCFMOB_ = std::unique_ptr<H5::DataSet>(new H5::DataSet(
        this->SCFGroup_->createDataSet("MOB",H5PredType<T>(),bdsp)));
    }
  }
};

template <typename T>
void SingleSlater<T>::fockDamping() {
  if(this->dampParam <= 0.) return;
  // Scalar
  //prettyPrintSmart(this->fileio_->out,*this->fockOrthoScalar_,"FockScalar predamp");
  this->SCFOrthoFScalar_->read(this->NBSqScratch_->data(),H5PredType<T>());
  *this->fockOrthoScalar_  *= (1-dampParam);
  //prettyPrintSmart(this->fileio_->out,*this->fockOrthoScalar_,"FockScalar scaled");
  *this->NBSqScratch_ *= dampParam;
  *this->fockOrthoScalar_  += *this->NBSqScratch_;
  //prettyPrintSmart(this->fileio_->out,*this->fockOrthoScalar_,"FockScalar damped");

  // Mz
  if(this->nTCS_ == 2 or !this->isClosedShell) {
    this->SCFOrthoFMz_->read(this->NBSqScratch_->data(),H5PredType<T>());
    *this->fockOrthoMz_  *= (1-dampParam);
    *this->NBSqScratch_ *= dampParam;
    *this->fockOrthoMz_  += *this->NBSqScratch_;

  }
  
  if(this->nTCS_ == 2){
    // My
    this->SCFOrthoFMy_->read(this->NBSqScratch_->data(),H5PredType<T>());
    *this->fockOrthoMy_  *= (1-dampParam);
    *this->NBSqScratch_ *= dampParam;
    *this->fockOrthoMy_  += *this->NBSqScratch_;

    // Mx
    this->SCFOrthoFMx_->read(this->NBSqScratch_->data(),H5PredType<T>());
    *this->fockOrthoMx_  *= (1-dampParam);
    *this->NBSqScratch_ *= dampParam;
    *this->fockOrthoMx_  += *this->NBSqScratch_;

  }

}

template <typename T>
void SingleSlater<T>::writeSCFFiles() {
  // Scalar
  this->SCFDensityScalar_->write(this->onePDMScalar_->data(),H5PredType<T>());
  this->SCFOrthoDScalar_->write(
    this->onePDMOrthoScalar_->data(),H5PredType<T>());

  this->SCFFockScalar_->write(this->fockScalar_->data(),H5PredType<T>());
  this->SCFOrthoFScalar_->write(this->fockOrthoScalar_->data(),H5PredType<T>());

  this->SCFPTScalar_->write(this->PTScalar_->data(),H5PredType<T>());
  this->SCFMOA_->write(this->moA_->data(),H5PredType<T>());

  if(this->nTCS_ == 2 or !this->isClosedShell) {
    // Mz
    this->SCFDensityMz_->write(this->onePDMMz_->data(),H5PredType<T>());
    this->SCFOrthoDMz_->write(this->onePDMOrthoMz_->data(),H5PredType<T>());

    this->SCFFockMz_->write(this->fockMz_->data(),H5PredType<T>());
    this->SCFOrthoFMz_->write(this->fockOrthoMz_->data(),H5PredType<T>());

    this->SCFPTMz_->write(this->PTMz_->data(),H5PredType<T>());
  }

  if(this->nTCS_ == 1 and !this->isClosedShell)
    this->SCFMOB_->write(this->moB_->data(),H5PredType<T>());

  if(this->nTCS_ == 2) {
    // Mx
    this->SCFDensityMx_->write(this->onePDMMx_->data(),H5PredType<T>());
    this->SCFOrthoDMx_->write(this->onePDMOrthoMx_->data(),H5PredType<T>());

    this->SCFFockMx_->write(this->fockMx_->data(),H5PredType<T>());
    this->SCFOrthoFMx_->write(this->fockOrthoMx_->data(),H5PredType<T>());

    this->SCFPTMx_->write(this->PTMx_->data(),H5PredType<T>());



    // My
    this->SCFDensityMy_->write(this->onePDMMy_->data(),H5PredType<T>());
    this->SCFOrthoDMy_->write(this->onePDMOrthoMy_->data(),H5PredType<T>());

    this->SCFFockMy_->write(this->fockMy_->data(),H5PredType<T>());
    this->SCFOrthoFMy_->write(this->fockOrthoMy_->data(),H5PredType<T>());

    this->SCFPTMy_->write(this->PTMy_->data(),H5PredType<T>());
  }
};
