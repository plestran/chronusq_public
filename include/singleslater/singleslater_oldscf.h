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

template<typename T>
void SingleSlater<T>::initSCFMem2(){
  this->initSCFPtr();


  auto NTCSxNBASIS = this->nTCS_ * this->nBasis_;
  auto NSQ = NTCSxNBASIS  * NTCSxNBASIS;

//this->LWORK_     = 5 * std::max(NTCSxNBASIS,this->nDIISExtrap_);
//this->LRWORK_    = 3 * std::max(NTCSxNBASIS,this->nDIISExtrap_);
  this->lenCoeff_  = this->nDIISExtrap_;

  
  this->POldAlphaMem_  = this->memManager_->template malloc<T>(NSQ);
  this->ErrorAlphaMem_ = 
    this->memManager_->template malloc<T>(NSQ*(this->nDIISExtrap_ -1));
  this->FADIIS_        = 
    this->memManager_->template malloc<T>(NSQ*(this->nDIISExtrap_ -1));

  if(this->nTCS_ == 1 && !this->isClosedShell) {
    this->POldBetaMem_  = this->memManager_->template malloc<T>(NSQ);
    this->ErrorBetaMem_ = 
      this->memManager_->template malloc<T>(NSQ*(this->nDIISExtrap_ -1));
    this->FBDIIS_       = 
      this->memManager_->template malloc<T>(NSQ*(this->nDIISExtrap_ -1));
  }

  if(this->Ref_ == CUHF) {
    this->delFMem_   = this->memManager_->template malloc<T>(NSQ);
    this->lambdaMem_ = this->memManager_->template malloc<T>(NSQ);
    this->PNOMem_    = this->memManager_->template malloc<T>(NSQ);
    this->occNumMem_ = this->memManager_->template malloc<double>(NTCSxNBASIS);
  }

//this->WORK_  = this->memManager_->template malloc<T>(this->LWORK_);
//if(typeid(T).hash_code() == typeid(dcomplex).hash_code())
//  this->RWORK_ = this->memManager_->template malloc<double>(this->LRWORK_);
  
  // New DIIS Files
  if(!this->doDIIS) return;

  std::vector<hsize_t> dims;
  dims.push_back(this->nDIISExtrap_ - 1);
  dims.push_back(this->nBasis_);
  dims.push_back(this->nBasis_);

  this->FScalarDIIS_ = 
    this->fileio_->createScratchPartition(H5PredType<T>(),
      "Fock (Ortho_Scalar) For DIIS Extrapoloation",dims);
  this->DScalarDIIS_ = 
    this->fileio_->createScratchPartition(H5PredType<T>(),
      "Density (Ortho_Scalar) For DIIS Extrapoloation",dims);
  this->EScalarDIIS_ = 
    this->fileio_->createScratchPartition(H5PredType<T>(),
      "Error Metric [F,D] (Scalar) For DIIS Extrapoloation",dims);

  if(this->nTCS_ == 2 || !this->isClosedShell) {
    this->FMzDIIS_ = 
      this->fileio_->createScratchPartition(H5PredType<T>(),
        "Fock (Ortho_Mz) For DIIS Extrapoloation",dims);
    this->DMzDIIS_ = 
      this->fileio_->createScratchPartition(H5PredType<T>(),
        "Density (Ortho_Mz) For DIIS Extrapoloation",dims);
    this->EMzDIIS_ = 
      this->fileio_->createScratchPartition(H5PredType<T>(),
        "Error Metric [F,D] (Mz) For DIIS Extrapoloation",dims);
  }

  if(this->nTCS_ == 2) {
    this->FMyDIIS_ = 
      this->fileio_->createScratchPartition(H5PredType<T>(),
        "Fock (Ortho_My) For DIIS Extrapoloation",dims);
    this->DMyDIIS_ = 
      this->fileio_->createScratchPartition(H5PredType<T>(),
        "Density (Ortho_My) For DIIS Extrapoloation",dims);
    this->EMyDIIS_ = 
      this->fileio_->createScratchPartition(H5PredType<T>(),
        "Error Metric [F,D] (My) For DIIS Extrapoloation",dims);

    this->FMxDIIS_ = 
      this->fileio_->createScratchPartition(H5PredType<T>(),
        "Fock (Ortho_Mx) For DIIS Extrapoloation",dims);
    this->DMxDIIS_ = 
      this->fileio_->createScratchPartition(H5PredType<T>(),
        "Density (Ortho_Mx) For DIIS Extrapoloation",dims);
    this->EMxDIIS_ = 
      this->fileio_->createScratchPartition(H5PredType<T>(),
        "Error Metric [F,D] (Mx) For DIIS Extrapoloation",dims);
  }

}; //initSCFMem2


template<typename T>
void SingleSlater<T>::SCF2(){
  // Compute the 1-Body Integrals if they haven't already been computed
  if(!this->aointegrals_->haveAOOneE) this->aointegrals_->computeAOOneE();
  if(this->printLevel_ > 0)
    this->printSCFHeader(this->fileio_->out);

  // Allocate SCF Memory
  this->initSCFMem2();

  bool doLevelShift; // dynamic boolean to decide when to toggle level shifting
  size_t iter;

  // FIXME: DIIS for CUHF NYI
  if(this->Ref_ == CUHF) this->doDIIS = false;
  
/*
  // Populate Orthonormal densities
  // Po(0)
  this->orthoDen2();
*/

  this->fileio_->out << "    *** INITIAL GUESS ENERGY = " << this->totalEnergy << " Eh ***" << endl;
  // SCF Iterations
  for(iter = 0; iter < this->maxSCFIter_; iter++){
    auto SCFStart = std::chrono::high_resolution_clock::now();

    // Copy the density to allow evaluation of RMS 
    // POLD = P(K)
    this->copyDen();

    // If this is CUHF (Scuseria), form the NOs and retransform the Fock matrix
    // accordingly
    if(this->Ref_ == CUHF) {
      this->formNO();
      this->fockCUHF();
    }

    // Transform (all of) the Fock matri(x,cies) from the AO to the 
    // orthonormal basis using the transformation matrix stored in AOIntegrals
    // Fo(K)
    this->orthoFock();

    // Copy the orthonormal Fock matrix to disk for DIIS
    if(this->doDIIS && iter >= this->iDIISStart_){
      if(iter == this->iDIISStart_ && this->printLevel_ > 0)
        this->fileio_->out << 
          std::setw(2) << " " <<
          std::setw(4) << " " <<
          "*** Starting DIIS ***" << endl;
      // Write Fo(K) to disk
      this->cpyOrthoFock2(iter - this->iDIISStart_); 
      // Write Po(K) to disk
      this->cpyOrthoDen2(iter - this->iDIISStart_); 
      // Write E(K) = [Fo(K),Po(K)] to disk
      this->genDIISCom(iter - this->iDIISStart_);
/*
      // Extrapolate Fo(K) -> C(I)Fo(I) using CDIIS
      if((iter - this->iDIISStart_) % (this->nDIISExtrap_-1) == 
         (this->nDIISExtrap_-2) && iter != 0) 
        this->CDIIS2();
*/
    }

//JJGS
    if(this->doITP) {
      if(iter == 0) this->diagFock2(); // Need one diagonalization to init guess
      this->doImagTimeProp(this->dt);
    } else {
      doLevelShift = iter >= this->iStartLevelShift_ ;
      doLevelShift = doLevelShift && 
        (iter < this->iStartLevelShift_ + this->nLevelShift_);
  
      // Fo(K) = Fo(K) + L
      if(doLevelShift) this->levelShift2();  //Level-shift for AO basis
      // Obtains Co(K+1)
      this->diagFock2();
    }
//JJGE

    // Optionally mix the orbitals (FIXME: This needs to be addressed)
    if(iter == 0 && this->guess_ != READ) this->mixOrbitalsSCF();


    // Form the density(s) in the orthonormal basis
    // P(K+1) = Co(K+1)*Co(K+1)**H
    this->formDensity();
  
    // Copy the data from onePDM? -> onePDMOrtho? (calculated above)
    // Po(K+1) = P(K+1)
    this->cpyAOtoOrthoDen();

    // Transform the density(s) back into the AO basis
    //this->orthoDen();
    // Po(K+1) -> P(K+1)
    this->unOrthoDen();

    // Form the AO Fock matrix(s)
    // F(K+1) = H + G[P(K+1)]
    this->formFock();

    // Check if interrupt has been encountered from python API
    if(PyErr_CheckSignals() == -1)
      CErr("Keyboard Interrupt in SCF!",this->fileio_->out);

/*
    // Perform DIIS Extrapolation for iterations greater than iDiisStart_
    //   This if-block both copies the Fock and error vectors  into scratch space
    //   for use with CDIIS as well as performing the extrapolation every nDIISExtrap
    //   steps
    //
    if(this->doDIIS && iter >= this->iDIISStart_){ 

      // Compute the error metric [FDS,SPF]
      this->genDComm2(iter - this->iDIISStart_);
      // Copy the Fock matrix(s) over for the extrapolation
      //   FIXME: Is this correct?? 
      //   Nakamo states that we should be extrapolating the Density not Fock
      this->CpyFock(iter - this->iDIISStart_);   

      // Perform the DIIS extrapolation every nDIISExtrap steps
      if((iter - this->iDIISStart_) % (this->nDIISExtrap_-1) == 
         (this->nDIISExtrap_-2) && iter != 0) 
//        this->CDIIS();
        this->CDIIS2();
    }
*/

    // Evaluate convergence critera
    this->evalConver(iter);
    this->nSCFIter++;

    auto SCFEnd = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> SCFD = SCFEnd - SCFStart;
//  cout << "SCF Time (" << iter << ") = " << SCFD.count() << endl; 
    
    // Break if converged
    if(this->isConverged) break;
  };
  // WARNING: MO Coefficients are not transformed to and from the
  // orthonormal basis throughout the SCF and must be transformed
  // back at the end for and post SCF to be functional
  this->backTransformMOs();

  this->cleanupSCFMem2();
  this->fixPhase();

  if(!this->isConverged)
    CErr("SCF Failed to converge within maximum number of iterations",
        this->fileio_->out);

  if(this->printLevel_ > 0 ){
    if(this->isConverged){
      this->fileio_->out 
        << endl << "SCF Completed: E(" << this->SCFTypeShort_ << ") = ";
      this->fileio_->out 
        << std::fixed << std::setprecision(10) 
        << this->totalEnergy << "  Eh after  " << iter + 1 
        << "  SCF Iterations" << endl;
    }
    this->fileio_->out << bannerEnd <<endl;
  }

};

template<typename T>
void SingleSlater<T>::cleanupSCFMem2(){
                                  
  auto NTCSxNBASIS = this->nTCS_ * this->nBasis_;
  auto NSQ = NTCSxNBASIS  * NTCSxNBASIS;

  this->memManager_->free(this->POldAlphaMem_,NSQ);  
  this->memManager_->free(this->ErrorAlphaMem_,(this->nDIISExtrap_-1)*NSQ); 
  this->memManager_->free(this->FADIIS_,(this->nDIISExtrap_-1)*NSQ);        
  if(this->nTCS_ == 1 && !this->isClosedShell) { 
    this->memManager_->free(this->POldBetaMem_,NSQ);  
    this->memManager_->free(this->ErrorBetaMem_,(this->nDIISExtrap_-1)*NSQ); 
    this->memManager_->free(this->FBDIIS_,(this->nDIISExtrap_-1)*NSQ);       
  }

  if(this->Ref_ == CUHF){ 
    this->memManager_->free(this->delFMem_,NSQ);   
    this->memManager_->free(this->lambdaMem_,NSQ); 
    this->memManager_->free(this->PNOMem_,NSQ);    
    this->memManager_->free(this->occNumMem_,NTCSxNBASIS); 
  }

//this->memManager_->free(this->WORK_,this->LWORK_);  
//if(typeid(T).hash_code() == typeid(dcomplex).hash_code()) 
//  this->memManager_->free(this->RWORK_,this->LRWORK_); 
  
  
}; //cleanupSCFMem2
