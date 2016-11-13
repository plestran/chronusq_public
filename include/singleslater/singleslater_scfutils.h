/* 
 *  This file is part of the Chronus Quantum (ChronusQ) software package
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
  this->initSCFFiles();

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

  hsize_t dims[2] = {this->nBasis_,this->nBasis_};
  hsize_t offset[2] = {0,0};

  H5::DataSpace memSpace(2,dims,NULL);
  H5::DataSpace curSpace = this->SCFDensityScalar_->getSpace();
  curSpace.selectHyperslab(H5S_SELECT_SET,dims,offset);

  // Scalar density RMS difference
  DeltaDScalar_->read(this->NBSqScratch_->data(),H5PredType<T>(),memSpace,
    curSpace);
  PSRMS = this->NBSqScratch_->norm();

  // Magnetization density RMS
     
  // || Pz(M) - Pz(M-1) ||^2
  if(this->nTCS_ == 2 or !this->isClosedShell){
    DeltaDMz_->read(this->NBSqScratch_->data(),H5PredType<T>(),memSpace,
      curSpace);
    PMRMS = this->NBSqScratch_->squaredNorm();
  }

  // || Py(M) - Py(M-1) ||^2 + || Px(M) - Px(M-1) ||^2
  if(this->nTCS_ == 2) {
    DeltaDMy_->read(this->NBSqScratch_->data(),H5PredType<T>(),memSpace,
      curSpace);
    PMRMS += this->NBSqScratch_->squaredNorm();
    DeltaDMx_->read(this->NBSqScratch_->data(),H5PredType<T>(),memSpace,
      curSpace);
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
  hsize_t dims[2] = {this->nBasis_,this->nBasis_};
  hsize_t offset[2] = {0,0};

  H5::DataSpace memSpace(2,dims,NULL);
  H5::DataSpace curSpace = this->SCFDensityScalar_->getSpace();
  curSpace.selectHyperslab(H5S_SELECT_SET,dims,offset);

  SCFDensityScalar_->read(this->NBSqScratch_->data(),H5PredType<T>(),
    memSpace,curSpace);
  (*this->NBSqScratch2_) = (*this->onePDMScalar_) - (*this->NBSqScratch_);
  DeltaDScalar_->write(this->NBSqScratch2_->data(),H5PredType<T>(),
    memSpace,curSpace);

  if(this->nTCS_ == 2 or !this->isClosedShell){
    SCFDensityMz_->read(this->NBSqScratch_->data(),H5PredType<T>(),
      memSpace,curSpace);
    (*this->NBSqScratch2_) = (*this->onePDMMz_) - (*this->NBSqScratch_);
    DeltaDMz_->write(this->NBSqScratch2_->data(),H5PredType<T>(),
      memSpace,curSpace);
  }

  if(this->nTCS_ == 2) {
    SCFDensityMy_->read(this->NBSqScratch_->data(),H5PredType<T>(),
      memSpace,curSpace);
    (*this->NBSqScratch2_) = (*this->onePDMMy_) - (*this->NBSqScratch_);
    DeltaDMy_->write(this->NBSqScratch2_->data(),H5PredType<T>(),
      memSpace,curSpace);

    SCFDensityMx_->read(this->NBSqScratch_->data(),H5PredType<T>(),
      memSpace,curSpace);
    (*this->NBSqScratch2_) = (*this->onePDMMx_) - (*this->NBSqScratch_);
    DeltaDMx_->write(this->NBSqScratch2_->data(),H5PredType<T>(),
      memSpace,curSpace);
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
  hsize_t dims[2] = {this->nBasis_,this->nBasis_};
  hsize_t offset[2] = {0,0};

  H5::DataSpace memSpace(2,dims,NULL);
  H5::DataSpace curSpace = this->SCFDensityScalar_->getSpace();
  curSpace.selectHyperslab(H5S_SELECT_SET,dims,offset);

  DeltaDScalar_->read(this->onePDMScalar_->data(),H5PredType<T>(),
    memSpace,curSpace);
  if(this->nTCS_ == 2 or !this->isClosedShell){
    DeltaDMz_->read(this->onePDMMz_->data(),H5PredType<T>(),memSpace,curSpace);
  }
  if(this->nTCS_ == 2) {
    DeltaDMy_->read(this->onePDMMy_->data(),H5PredType<T>(),memSpace,curSpace);
    DeltaDMx_->read(this->onePDMMx_->data(),H5PredType<T>(),memSpace,curSpace);
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
  hsize_t dims[2] = {this->nBasis_,this->nBasis_};
  hsize_t offset[2] = {0,0};

  H5::DataSpace memSpace(2,dims,NULL);
  H5::DataSpace curSpace = this->SCFDensityScalar_->getSpace();
  curSpace.selectHyperslab(H5S_SELECT_SET,dims,offset);

  SCFDensityScalar_->read(this->onePDMScalar_->data(),H5PredType<T>(),
    memSpace,curSpace);
  if(this->nTCS_ == 2 or !this->isClosedShell){
    SCFDensityMz_->read(this->onePDMMz_->data(),H5PredType<T>(),memSpace,
      curSpace);
  }
  if(this->nTCS_ == 2) {
    SCFDensityMy_->read(this->onePDMMy_->data(),H5PredType<T>(),memSpace,
      curSpace);
    SCFDensityMx_->read(this->onePDMMx_->data(),H5PredType<T>(),memSpace,
      curSpace);
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

//if(!this->isPrimary) return;
  std::vector<hsize_t> dims;
  dims.push_back(this->basisset_->nBasis());
  dims.push_back(this->basisset_->nBasis());

  std::vector<hsize_t> bigDims;
  bigDims.push_back(this->nTCS_*this->basisset_->nBasis());
  bigDims.push_back(this->nTCS_*this->basisset_->nBasis());

  H5::DataSpace dsp(2,&dims[0]);
  H5::DataSpace bdsp(2,&bigDims[0]);

  // Stuff for extendable files
  hsize_t chunk_dims[2] = {2,2};
  H5::DSetCreatPropList prop;
  if(!this->isPrimary) prop.setChunk(2,chunk_dims);

  // Check if SCF group exists
  if(this->isPrimary) {
    try {
      this->SCFGroup_ = std::unique_ptr<H5::Group>(new H5::Group(
        this->fileio_->restart->openGroup("/SCF")));
    } catch(...) {
      this->SCFGroup_ = std::unique_ptr<H5::Group>(new H5::Group(
        this->fileio_->restart->createGroup("/SCF")));
    }
  } else {
    try {
      this->SCFGroup_ = std::unique_ptr<H5::Group>(new H5::Group(
        this->fileio_->scr->openGroup("/NONPRIMARY_SCF")));
    } catch(...) {
      this->SCFGroup_ = std::unique_ptr<H5::Group>(new H5::Group(
        this->fileio_->scr->createGroup("/NONPRIMARY_SCF")));
    }
  }


  // Density Files

  // Scalar
  try {
    this->SCFDensityScalar_ = std::unique_ptr<H5::DataSet>( new H5::DataSet(
      this->SCFGroup_->openDataSet("DensityScalar")));

    // Check if resize needed
    if(!this->isPrimary){
      std::vector<hsize_t> prev_dims(2);

      this->SCFDensityScalar_->getSpace().getSimpleExtentDims(&prev_dims[0]);
      if(prev_dims[0] < dims[0])
        this->SCFDensityScalar_->extend(&dims[0]);
    }
  } catch(...) {
    if(this->guess_ == READ) 
      CErr("FATAL: READ Guess cannot find DensityScalar",this->fileio_->out);

    this->SCFDensityScalar_ = std::unique_ptr<H5::DataSet>( new H5::DataSet(
      this->SCFGroup_->createDataSet("DensityScalar",H5PredType<T>(),dsp,prop)));
  }

  // Mz
  if(this->nTCS_ == 2 or !this->isClosedShell) {
    try {
      this->SCFDensityMz_ = std::unique_ptr<H5::DataSet>( new H5::DataSet(
        this->SCFGroup_->openDataSet("DensityMz")));

      // Check if resize needed
      if(!this->isPrimary){
        std::vector<hsize_t> prev_dims(2);
     
        this->SCFDensityMz_->getSpace().getSimpleExtentDims(&prev_dims[0]);
        if(prev_dims[0] < dims[0])
          this->SCFDensityMz_->extend(&dims[0]);
      }
    } catch(...) {
      this->SCFDensityMz_ = std::unique_ptr<H5::DataSet>( new H5::DataSet(
        this->SCFGroup_->createDataSet("DensityMz",H5PredType<T>(),dsp,prop)));
    }
  }

  // Mx / My
  if(this->nTCS_ == 2) {
    // Mx
    try {
      this->SCFDensityMx_ = std::unique_ptr<H5::DataSet>( new H5::DataSet(
        this->SCFGroup_->openDataSet("DensityMx")));

      // Check if resize needed
      if(!this->isPrimary){
        std::vector<hsize_t> prev_dims(2);
     
        this->SCFDensityMx_->getSpace().getSimpleExtentDims(&prev_dims[0]);
        if(prev_dims[0] < dims[0])
          this->SCFDensityMx_->extend(&dims[0]);
      }
    } catch(...) {
      this->SCFDensityMx_ = std::unique_ptr<H5::DataSet>( new H5::DataSet(
        this->SCFGroup_->createDataSet("DensityMx",H5PredType<T>(),dsp,prop)));
    }

    // My
    try {
      this->SCFDensityMy_ = std::unique_ptr<H5::DataSet>( new H5::DataSet(
        this->SCFGroup_->openDataSet("DensityMy")));

      // Check if resize needed
      if(!this->isPrimary){
        std::vector<hsize_t> prev_dims(2);
     
        this->SCFDensityMy_->getSpace().getSimpleExtentDims(&prev_dims[0]);
        if(prev_dims[0] < dims[0])
          this->SCFDensityMy_->extend(&dims[0]);
      }
    } catch(...) {
      this->SCFDensityMy_ = std::unique_ptr<H5::DataSet>( new H5::DataSet(
        this->SCFGroup_->createDataSet("DensityMy",H5PredType<T>(),dsp,prop)));
    }
  }

  // Ortho Density Files

  // Scalar
  try {
    this->SCFOrthoDScalar_ = std::unique_ptr<H5::DataSet>( new H5::DataSet(
      this->SCFGroup_->openDataSet("OrthoDScalar")));

    // Check if resize needed
    if(!this->isPrimary){
      std::vector<hsize_t> prev_dims(2);
    
      this->SCFOrthoDScalar_->getSpace().getSimpleExtentDims(&prev_dims[0]);
      if(prev_dims[0] < dims[0])
        this->SCFOrthoDScalar_->extend(&dims[0]);
    }
  } catch(...) {
    this->SCFOrthoDScalar_ = std::unique_ptr<H5::DataSet>( new H5::DataSet(
      this->SCFGroup_->createDataSet("OrthoDScalar",H5PredType<T>(),dsp,prop)));
  }

  // Mz
  if(this->nTCS_ == 2 or !this->isClosedShell) {
    try {
      this->SCFOrthoDMz_ = std::unique_ptr<H5::DataSet>( new H5::DataSet(
        this->SCFGroup_->openDataSet("OrthoDMz")));

      // Check if resize needed
      if(!this->isPrimary){
        std::vector<hsize_t> prev_dims(2);
      
        this->SCFOrthoDMz_->getSpace().getSimpleExtentDims(&prev_dims[0]);
        if(prev_dims[0] < dims[0])
          this->SCFOrthoDMz_->extend(&dims[0]);
      }
    } catch(...) {
      this->SCFOrthoDMz_ = std::unique_ptr<H5::DataSet>( new H5::DataSet(
        this->SCFGroup_->createDataSet("OrthoDMz",H5PredType<T>(),dsp,prop)));
    }
  }

  // Mx / My
  if(this->nTCS_ == 2) {
    // Mx
    try {
      this->SCFOrthoDMx_ = std::unique_ptr<H5::DataSet>( new H5::DataSet(
        this->SCFGroup_->openDataSet("OrthoDMx")));

      // Check if resize needed
      if(!this->isPrimary){
        std::vector<hsize_t> prev_dims(2);
      
        this->SCFOrthoDMx_->getSpace().getSimpleExtentDims(&prev_dims[0]);
        if(prev_dims[0] < dims[0])
          this->SCFOrthoDMx_->extend(&dims[0]);
      }
    } catch(...) {
      this->SCFOrthoDMx_ = std::unique_ptr<H5::DataSet>( new H5::DataSet(
        this->SCFGroup_->createDataSet("OrthoDMx",H5PredType<T>(),dsp,prop)));
    }

    // My
    try {
      this->SCFOrthoDMy_ = std::unique_ptr<H5::DataSet>( new H5::DataSet(
        this->SCFGroup_->openDataSet("OrthoDMy")));

      // Check if resize needed
      if(!this->isPrimary){
        std::vector<hsize_t> prev_dims(2);
      
        this->SCFOrthoDMy_->getSpace().getSimpleExtentDims(&prev_dims[0]);
        if(prev_dims[0] < dims[0])
          this->SCFOrthoDMy_->extend(&dims[0]);
      }
    } catch(...) {
      this->SCFOrthoDMy_ = std::unique_ptr<H5::DataSet>( new H5::DataSet(
        this->SCFGroup_->createDataSet("OrthoDMy",H5PredType<T>(),dsp,prop)));
    }
  }



  // Fock Files

  // Scalar
  try {
    this->SCFFockScalar_ = std::unique_ptr<H5::DataSet>( new H5::DataSet(
      this->SCFGroup_->openDataSet("FockScalar")));

    // Check if resize needed
    if(!this->isPrimary){
      std::vector<hsize_t> prev_dims(2);

      this->SCFFockScalar_->getSpace().getSimpleExtentDims(&prev_dims[0]);
      if(prev_dims[0] < dims[0])
        this->SCFFockScalar_->extend(&dims[0]);
    }
  } catch(...) {
    this->SCFFockScalar_ = std::unique_ptr<H5::DataSet>( new H5::DataSet(
      this->SCFGroup_->createDataSet("FockScalar",H5PredType<T>(),dsp,prop)));
  }

  // Mz
  if(this->nTCS_ == 2 or !this->isClosedShell) {
    try {
      this->SCFFockMz_ = std::unique_ptr<H5::DataSet>( new H5::DataSet(
        this->SCFGroup_->openDataSet("FockMz")));

      // Check if resize needed
      if(!this->isPrimary){
        std::vector<hsize_t> prev_dims(2);
     
        this->SCFFockMz_->getSpace().getSimpleExtentDims(&prev_dims[0]);
        if(prev_dims[0] < dims[0])
          this->SCFFockMz_->extend(&dims[0]);
      }
    } catch(...) {
      this->SCFFockMz_ = std::unique_ptr<H5::DataSet>( new H5::DataSet(
        this->SCFGroup_->createDataSet("FockMz",H5PredType<T>(),dsp,prop)));
    }
  }

  // Mx / My
  if(this->nTCS_ == 2) {
    // Mx
    try {
      this->SCFFockMx_ = std::unique_ptr<H5::DataSet>( new H5::DataSet(
        this->SCFGroup_->openDataSet("FockMx")));

      // Check if resize needed
      if(!this->isPrimary){
        std::vector<hsize_t> prev_dims(2);
     
        this->SCFFockMx_->getSpace().getSimpleExtentDims(&prev_dims[0]);
        if(prev_dims[0] < dims[0])
          this->SCFFockMx_->extend(&dims[0]);
      }
    } catch(...) {
      this->SCFFockMx_ = std::unique_ptr<H5::DataSet>( new H5::DataSet(
        this->SCFGroup_->createDataSet("FockMx",H5PredType<T>(),dsp,prop)));
    }

    // My
    try {
      this->SCFFockMy_ = std::unique_ptr<H5::DataSet>( new H5::DataSet(
        this->SCFGroup_->openDataSet("FockMy")));

      // Check if resize needed
      if(!this->isPrimary){
        std::vector<hsize_t> prev_dims(2);
     
        this->SCFFockMy_->getSpace().getSimpleExtentDims(&prev_dims[0]);
        if(prev_dims[0] < dims[0])
          this->SCFFockMy_->extend(&dims[0]);
      }
    } catch(...) {
      this->SCFFockMy_ = std::unique_ptr<H5::DataSet>( new H5::DataSet(
        this->SCFGroup_->createDataSet("FockMy",H5PredType<T>(),dsp,prop)));
    }
  }

  // Ortho Fock Files

  // Scalar
  try {
    this->SCFOrthoFScalar_ = std::unique_ptr<H5::DataSet>( new H5::DataSet(
      this->SCFGroup_->openDataSet("OrthoFScalar")));

    // Check if resize needed
    if(!this->isPrimary){
      std::vector<hsize_t> prev_dims(2);

      this->SCFOrthoFScalar_->getSpace().getSimpleExtentDims(&prev_dims[0]);
      if(prev_dims[0] < dims[0])
        this->SCFOrthoFScalar_->extend(&dims[0]);
    }
  } catch(...) {
    this->SCFOrthoFScalar_ = std::unique_ptr<H5::DataSet>( new H5::DataSet(
      this->SCFGroup_->createDataSet("OrthoFScalar",H5PredType<T>(),dsp,prop)));
  }

  // Mz
  if(this->nTCS_ == 2 or !this->isClosedShell) {
    try {
      this->SCFOrthoFMz_ = std::unique_ptr<H5::DataSet>( new H5::DataSet(
        this->SCFGroup_->openDataSet("OrthoFMz")));

      // Check if resize needed
      if(!this->isPrimary){
        std::vector<hsize_t> prev_dims(2);
     
        this->SCFOrthoFMz_->getSpace().getSimpleExtentDims(&prev_dims[0]);
        if(prev_dims[0] < dims[0])
          this->SCFOrthoFMz_->extend(&dims[0]);
      }
    } catch(...) {
      this->SCFOrthoFMz_ = std::unique_ptr<H5::DataSet>( new H5::DataSet(
        this->SCFGroup_->createDataSet("OrthoFMz",H5PredType<T>(),dsp,prop)));
    }
  }

  // Mx / My
  if(this->nTCS_ == 2) {
    // Mx
    try {
      this->SCFOrthoFMx_ = std::unique_ptr<H5::DataSet>( new H5::DataSet(
        this->SCFGroup_->openDataSet("OrthoFMx")));

      // Check if resize needed
      if(!this->isPrimary){
        std::vector<hsize_t> prev_dims(2);
     
        this->SCFOrthoFMx_->getSpace().getSimpleExtentDims(&prev_dims[0]);
        if(prev_dims[0] < dims[0])
          this->SCFOrthoFMx_->extend(&dims[0]);
      }
    } catch(...) {
      this->SCFOrthoFMx_ = std::unique_ptr<H5::DataSet>( new H5::DataSet(
        this->SCFGroup_->createDataSet("OrthoFMx",H5PredType<T>(),dsp,prop)));
    }

    // My
    try {
      this->SCFOrthoFMy_ = std::unique_ptr<H5::DataSet>( new H5::DataSet(
        this->SCFGroup_->openDataSet("OrthoFMy")));

      // Check if resize needed
      if(!this->isPrimary){
        std::vector<hsize_t> prev_dims(2);
     
        this->SCFOrthoFMy_->getSpace().getSimpleExtentDims(&prev_dims[0]);
        if(prev_dims[0] < dims[0])
          this->SCFOrthoFMy_->extend(&dims[0]);
      }
    } catch(...) {
      this->SCFOrthoFMy_ = std::unique_ptr<H5::DataSet>( new H5::DataSet(
        this->SCFGroup_->createDataSet("OrthoFMy",H5PredType<T>(),dsp,prop)));
    }
  }

  // PT Files

  // Scalar
  try {
    this->SCFPTScalar_ = std::unique_ptr<H5::DataSet>( new H5::DataSet(
      this->SCFGroup_->openDataSet("PTScalar")));

    // Check if resize needed
    if(!this->isPrimary){
      std::vector<hsize_t> prev_dims(2);

      this->SCFPTScalar_->getSpace().getSimpleExtentDims(&prev_dims[0]);
      if(prev_dims[0] < dims[0])
        this->SCFPTScalar_->extend(&dims[0]);
    }
  } catch(...) {
    this->SCFPTScalar_ = std::unique_ptr<H5::DataSet>( new H5::DataSet(
      this->SCFGroup_->createDataSet("PTScalar",H5PredType<T>(),dsp,prop)));
  }

  // Mz
  if(this->nTCS_ == 2 or !this->isClosedShell) {
    try {
      this->SCFPTMz_ = std::unique_ptr<H5::DataSet>( new H5::DataSet(
        this->SCFGroup_->openDataSet("PTMz")));

      // Check if resize needed
      if(!this->isPrimary){
        std::vector<hsize_t> prev_dims(2);

        this->SCFPTMz_->getSpace().getSimpleExtentDims(&prev_dims[0]);
        if(prev_dims[0] < dims[0])
          this->SCFPTMz_->extend(&dims[0]);
      }
    } catch(...) {
      this->SCFPTMz_ = std::unique_ptr<H5::DataSet>( new H5::DataSet(
        this->SCFGroup_->createDataSet("PTMz",H5PredType<T>(),dsp,prop)));
    }
  }

  // Mx / My
  if(this->nTCS_ == 2) {
    // Mx
    try {
      this->SCFPTMx_ = std::unique_ptr<H5::DataSet>( new H5::DataSet(
        this->SCFGroup_->openDataSet("PTMx")));

      // Check if resize needed
      if(!this->isPrimary){
        std::vector<hsize_t> prev_dims(2);

        this->SCFPTMx_->getSpace().getSimpleExtentDims(&prev_dims[0]);
        if(prev_dims[0] < dims[0])
          this->SCFPTMx_->extend(&dims[0]);
      }
    } catch(...) {
      this->SCFPTMx_ = std::unique_ptr<H5::DataSet>( new H5::DataSet(
        this->SCFGroup_->createDataSet("PTMx",H5PredType<T>(),dsp,prop)));
    }

    // My
    try {
      this->SCFPTMy_ = std::unique_ptr<H5::DataSet>( new H5::DataSet(
        this->SCFGroup_->openDataSet("PTMy")));

      // Check if resize needed
      if(!this->isPrimary){
        std::vector<hsize_t> prev_dims(2);

        this->SCFPTMy_->getSpace().getSimpleExtentDims(&prev_dims[0]);
        if(prev_dims[0] < dims[0])
          this->SCFPTMy_->extend(&dims[0]);
      }
    } catch(...) {
      this->SCFPTMy_ = std::unique_ptr<H5::DataSet>( new H5::DataSet(
        this->SCFGroup_->createDataSet("PTMy",H5PredType<T>(),dsp,prop)));
    }
  }

  // MO Files
  
  // MO (TCS) / MOA
  try {
    this->SCFMOA_ = std::unique_ptr<H5::DataSet>(new H5::DataSet(
      this->SCFGroup_->openDataSet("MOA")));

    // Check if resize needed
    if(!this->isPrimary){
      std::vector<hsize_t> prev_dims(2);

      this->SCFMOA_->getSpace().getSimpleExtentDims(&prev_dims[0]);
      if(prev_dims[0] < bigDims[0])
        this->SCFMOA_->extend(&bigDims[0]);
    }
  } catch(...) {
    this->SCFMOA_ = std::unique_ptr<H5::DataSet>(new H5::DataSet(
      this->SCFGroup_->createDataSet("MOA",H5PredType<T>(),bdsp,prop)));
  }

  // MOB
  if(this->nTCS_ == 1 and !this->isClosedShell){
    try {
      this->SCFMOB_ = std::unique_ptr<H5::DataSet>(new H5::DataSet(
        this->SCFGroup_->openDataSet("MOB")));

      // Check if resize needed
      if(!this->isPrimary){
        std::vector<hsize_t> prev_dims(2);
     
        this->SCFMOB_->getSpace().getSimpleExtentDims(&prev_dims[0]);
        if(prev_dims[0] < bigDims[0])
          this->SCFMOB_->extend(&bigDims[0]);
      }
    } catch(...) {
      this->SCFMOB_ = std::unique_ptr<H5::DataSet>(new H5::DataSet(
        this->SCFGroup_->createDataSet("MOB",H5PredType<T>(),bdsp,prop)));
    }
  }
};

template <typename T>
void SingleSlater<T>::fockDamping() {
  if(this->dampParam <= 0.) return;
  hsize_t dims[2] = {this->nBasis_,this->nBasis_};
  hsize_t offset[2] = {0,0};

  H5::DataSpace memSpace(2,dims,NULL);
  H5::DataSpace curSpace = this->SCFDensityScalar_->getSpace();
  curSpace.selectHyperslab(H5S_SELECT_SET,dims,offset);

  // Scalar
  //prettyPrintSmart(this->fileio_->out,*this->fockOrthoScalar_,"FockScalar predamp");
  this->SCFOrthoFScalar_->read(this->NBSqScratch_->data(),H5PredType<T>(),
    memSpace,curSpace);
  *this->fockOrthoScalar_  *= (1-dampParam);
  //prettyPrintSmart(this->fileio_->out,*this->fockOrthoScalar_,"FockScalar scaled");
  *this->NBSqScratch_ *= dampParam;
  *this->fockOrthoScalar_  += *this->NBSqScratch_;
  //prettyPrintSmart(this->fileio_->out,*this->fockOrthoScalar_,"FockScalar damped");

  // Mz
  if(this->nTCS_ == 2 or !this->isClosedShell) {
    this->SCFOrthoFMz_->read(this->NBSqScratch_->data(),H5PredType<T>(),
      memSpace,curSpace);
    *this->fockOrthoMz_  *= (1-dampParam);
    *this->NBSqScratch_ *= dampParam;
    *this->fockOrthoMz_  += *this->NBSqScratch_;

  }
  
  if(this->nTCS_ == 2){
    // My
    this->SCFOrthoFMy_->read(this->NBSqScratch_->data(),H5PredType<T>(),
      memSpace,curSpace);
    *this->fockOrthoMy_  *= (1-dampParam);
    *this->NBSqScratch_ *= dampParam;
    *this->fockOrthoMy_  += *this->NBSqScratch_;

    // Mx
    this->SCFOrthoFMx_->read(this->NBSqScratch_->data(),H5PredType<T>(),
      memSpace,curSpace);
    *this->fockOrthoMx_  *= (1-dampParam);
    *this->NBSqScratch_ *= dampParam;
    *this->fockOrthoMx_  += *this->NBSqScratch_;

  }

}

template <typename T>
void SingleSlater<T>::writeSCFFiles() {
  hsize_t dims[2] = {this->nBasis_,this->nBasis_};
  hsize_t offset[2] = {0,0};

  H5::DataSpace memSpace(2,dims,NULL);
  H5::DataSpace curSpace = this->SCFDensityScalar_->getSpace();
  curSpace.selectHyperslab(H5S_SELECT_SET,dims,offset);

  // Scalar
  this->SCFDensityScalar_->write(this->onePDMScalar_->data(),H5PredType<T>(),
    memSpace,curSpace);
  this->SCFOrthoDScalar_->write(
    this->onePDMOrthoScalar_->data(),H5PredType<T>(),memSpace,curSpace);

  this->SCFFockScalar_->write(this->fockScalar_->data(),H5PredType<T>(),
    memSpace,curSpace);
  this->SCFOrthoFScalar_->write(this->fockOrthoScalar_->data(),H5PredType<T>(),
    memSpace,curSpace);

  this->SCFPTScalar_->write(this->PTScalar_->data(),H5PredType<T>(),
    memSpace,curSpace);
  this->SCFMOA_->write(this->moA_->data(),H5PredType<T>(),memSpace,curSpace);

  if(this->nTCS_ == 2 or !this->isClosedShell) {
    // Mz
    this->SCFDensityMz_->write(this->onePDMMz_->data(),H5PredType<T>(),
      memSpace,curSpace);
    this->SCFOrthoDMz_->write(this->onePDMOrthoMz_->data(),H5PredType<T>(),
      memSpace,curSpace);

    this->SCFFockMz_->write(this->fockMz_->data(),H5PredType<T>(),
      memSpace,curSpace);
    this->SCFOrthoFMz_->write(this->fockOrthoMz_->data(),H5PredType<T>(),
      memSpace,curSpace);

    this->SCFPTMz_->write(this->PTMz_->data(),H5PredType<T>(),
      memSpace,curSpace);
  }

  if(this->nTCS_ == 1 and !this->isClosedShell)
    this->SCFMOB_->write(this->moB_->data(),H5PredType<T>(),
      memSpace,curSpace);

  if(this->nTCS_ == 2) {
    // Mx
    this->SCFDensityMx_->write(this->onePDMMx_->data(),H5PredType<T>(),
      memSpace,curSpace);
    this->SCFOrthoDMx_->write(this->onePDMOrthoMx_->data(),H5PredType<T>(),
      memSpace,curSpace);

    this->SCFFockMx_->write(this->fockMx_->data(),H5PredType<T>(),
      memSpace,curSpace);
    this->SCFOrthoFMx_->write(this->fockOrthoMx_->data(),H5PredType<T>(),
      memSpace,curSpace);

    this->SCFPTMx_->write(this->PTMx_->data(),H5PredType<T>(),
      memSpace,curSpace);



    // My
    this->SCFDensityMy_->write(this->onePDMMy_->data(),H5PredType<T>(),
      memSpace,curSpace);
    this->SCFOrthoDMy_->write(this->onePDMOrthoMy_->data(),H5PredType<T>(),
      memSpace,curSpace);

    this->SCFFockMy_->write(this->fockMy_->data(),H5PredType<T>(),
      memSpace,curSpace);
    this->SCFOrthoFMy_->write(this->fockOrthoMy_->data(),H5PredType<T>(),
      memSpace,curSpace);

    this->SCFPTMy_->write(this->PTMy_->data(),H5PredType<T>(),
      memSpace,curSpace);
  }
};
