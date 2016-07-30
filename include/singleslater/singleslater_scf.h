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

template <typename T>
void SingleSlater<T>::initSCFPtr(){
  this->occNumMem_     = NULL;
//this->RWORK_         = NULL;

  this->FpAlphaMem_    = NULL;
  this->FpBetaMem_     = NULL;
  this->POldAlphaMem_  = NULL;
  this->POldBetaMem_   = NULL;
  this->ErrorAlphaMem_ = NULL;
  this->ErrorBetaMem_  = NULL;
  this->FADIIS_        = NULL;
  this->FBDIIS_        = NULL;
//this->WORK_          = NULL;
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
      "Fock (Scalar) For DIIS Extrapoloation",dims);
  this->DScalarDIIS_ = 
    this->fileio_->createScratchPartition(H5PredType<T>(),
      "Density (Scalar) For DIIS Extrapoloation",dims);
  this->DScalarDIIS_ = 
    this->fileio_->createScratchPartition(H5PredType<T>(),
      "Error Metric [F,D] (Scalar) For DIIS Extrapoloation",dims);

  if(this->nTCS_ == 2 || !this->isClosedShell) {
    this->FMzDIIS_ = 
      this->fileio_->createScratchPartition(H5PredType<T>(),
        "Fock (Mz) For DIIS Extrapoloation",dims);
    this->DMzDIIS_ = 
      this->fileio_->createScratchPartition(H5PredType<T>(),
        "Density (Mz) For DIIS Extrapoloation",dims);
    this->DMzDIIS_ = 
      this->fileio_->createScratchPartition(H5PredType<T>(),
        "Error Metric [F,D] (Mz) For DIIS Extrapoloation",dims);
  }

  if(this->nTCS_ == 2) {
    this->FMyDIIS_ = 
      this->fileio_->createScratchPartition(H5PredType<T>(),
        "Fock (My) For DIIS Extrapoloation",dims);
    this->DMyDIIS_ = 
      this->fileio_->createScratchPartition(H5PredType<T>(),
        "Density (My) For DIIS Extrapoloation",dims);
    this->DMyDIIS_ = 
      this->fileio_->createScratchPartition(H5PredType<T>(),
        "Error Metric [F,D] (My) For DIIS Extrapoloation",dims);

    this->FMxDIIS_ = 
      this->fileio_->createScratchPartition(H5PredType<T>(),
        "Fock (Mx) For DIIS Extrapoloation",dims);
    this->DMxDIIS_ = 
      this->fileio_->createScratchPartition(H5PredType<T>(),
        "Density (Mx) For DIIS Extrapoloation",dims);
    this->DMxDIIS_ = 
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

  bool doLevelShift; // dynamic boolean to decide when to turn off/on level shifting
  size_t iter;

  // SCF Iterations
  for(iter = 0; iter < this->maxSCFIter_; iter++){
    auto SCFStart = std::chrono::high_resolution_clock::now();

    // If this is CUHF (Scuseria), form the NOs and retransform the Fock matrix
    // accordingly
    if(this->Ref_ == CUHF) {
      this->formNO();
      this->fockCUHF();
    }

    // Transform (all of) the Fock matri(x,cies) from the AO to the orthonormal basis
    // using the transformation matrix stored in AOIntegrals
    this->orthoFock();

//JJGS
    if(this->doITP) {
      if(iter == 0) this->diagFock2(); // Need one diagonalization to init guess
      this->doImagTimeProp(this->dt);
    } else {
      doLevelShift = iter >= this->iStartLevelShift_ ;
      doLevelShift = doLevelShift && 
        (iter < this->iStartLevelShift_ + this->nLevelShift_);
  
      if(doLevelShift)
        this->levelShift2();	//Level-shift for AO basis
      this->diagFock2();
    }
//JJGE

    // Optionally mix the orbitals (FIXME: This needs to be addressed)
    if(iter == 0 && this->guess_ != READ) this->mixOrbitalsSCF();

    // Copy the density to allow evaluation of RMS (FIXME: this should just
    // copy the densities into scratch space as with the Focks to later evaluate
    // the DIIS criteria, this facilities EDIIS as well)
    this->copyDen();

    // Form the density(s) in the orthonormal basis
    this->formDensity();

    // Transform the density(s) back into the AO basis
    // FIXME: The name of this function needs to be changed, it is a misnomer.
    this->orthoDen();

    // Form the AO Fock matrix(s)
    this->formFock();

    // Check if interrupt has been encountered from python API
    if(PyErr_CheckSignals() == -1)
      CErr("Keyboard Interrupt in SCF!",this->fileio_->out);

    // Perform DIIS Extrapolation for iterations greater than iDiisStart_
    //   This if-block both copies the Fock and error vectors  into scratch space
    //   for use with CDIIS as well as performing the extrapolation every nDIISExtrap
    //   steps
    //
    // FIXME: DIIS for CUHF NYI
    if(this->Ref_ != CUHF && this->doDIIS && iter >= this->iDIISStart_){ 

      if(iter == this->iDIISStart_ && this->printLevel_ > 0)
        this->fileio_->out << 
          std::setw(2) << " " <<
          std::setw(4) << " " <<
          "*** Starting DIIS ***" << endl;
      // Compute the error metric [FDS,SPF]
      this->genDComm2(iter - this->iDIISStart_);
      // Copy the Fock matrix(s) over for the extrapolation
      //   FIXME: Is this correct?? 
      //   Nakamo states that we should be extrapolating the Density not Fock
      this->CpyFock(iter - this->iDIISStart_);   

      // Perform the DIIS extrapolation every nDIISExtrap steps
      if((iter - this->iDIISStart_) % (this->nDIISExtrap_-1) == 
         (this->nDIISExtrap_-2) && iter != 0) 
        this->CDIIS();
    }

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

template <typename T>
void SingleSlater<T>::copyDen(){
  auto NTCSxNBASIS = this->nTCS_*this->nBasis_;
  TMap POldAlpha(this->POldAlphaMem_,NTCSxNBASIS,NTCSxNBASIS);
  POldAlpha = (*this->onePDMA_);

  if(this->nTCS_ == 1 && !this->isClosedShell){
    TMap POldBeta(this->POldBetaMem_,NTCSxNBASIS,NTCSxNBASIS);
    POldBeta = (*this->onePDMB_);
  };
};
//DBWY and JJRADLER  -- Needs to be rewritten and have conditional for convergence added
//JJRADLER - This function should be called within SCF loop until convergence criterion is met.
//This algorithm follows Mehrotra
template<typename T>
void SingleSlater<T>::levelShift(){
//Construction of Delta matrix
// cout << "In Level Shift" << endl;
  double b = 2.42;	// 2 + .The meaning of the Universe.
  TMatrix deltaA = TMatrix::Zero(this->nBasis_, this->nBasis_);
  for(auto iVir = this->nOccA_; iVir < this->nBasis_; iVir++){
    deltaA(iVir, iVir) = b;
  }
  //Transformation of Delta matrix from MO to AO basis
  TMatrix transC = (*this->moA_).transpose();
  TMatrix dagC = transC.conjugate();
  TMatrix deltaAPrime = dagC * deltaA * (*this->moA_);	//Transformation of delta (MO) to deltaPrime (AO)
  
  //Shift current Fock matrix (F + Delta')
  for(auto iVir = this->nOccA_; iVir < this->nBasis_; iVir++){
    (*this->fockA_)(iVir, iVir) += deltaAPrime(iVir, iVir);	//Sums the diagonal elements of the vir-vir blocks only
  }

  //If CUHF is performed, apply same treatment to Beta MOs
  if(this->Ref_ == CUHF) {
    TMatrix deltaB = TMatrix::Zero(this->nBasis_, this->nBasis_);
    for(auto iVir = this->nOccB_; iVir < this->nBasis_; iVir++){
      deltaB(iVir, iVir) = b;
    }
    //Transformation of Delta matrix from MO to AO basis
    TMatrix transC = (*this->moB_).transpose();
    TMatrix dagC = transC.conjugate();
    TMatrix deltaBPrime = dagC * deltaB * (*this->moB_);	//Transformation of delta (MO) to deltaPrime (AO)
  
    //Shift current Fock matrix (F + Delta')
    for(auto iVir = this->nOccB_; iVir < this->nBasis_; iVir++){
      (*this->fockB_)(iVir, iVir) += deltaBPrime(iVir, iVir);	//Sums the diagonal elements of the vir-vir blocks only
    }
  }
}

//And this one is even simpler -- JJR
template<typename T>
void SingleSlater<T>::levelShift2(){
//double b = 2.42;	//2+.The meaning of The Universe

  T* FockA, *FockB;
  T* MOAVir, *MOBVir;
  T MOAMu, MOBMu;
//T* FockA = this->fockA_->data();
  FockA = this->fockOrthoA_->data();
  if(this->nTCS_ == 1 && !this->isClosedShell) 
    FockB = this->fockOrthoB_->data();
  
  for(auto iVir = this->nOccA_; iVir < this->nBasis_; iVir++){
    MOAVir = this->moA_->data() + iVir * this->nBasis_;
    if(this->nTCS_ == 1 && !this->isClosedShell) 
      MOBVir = this->moB_->data() + iVir * this->nBasis_;

    for(auto mu = 0; mu < this->nBasis_; mu++){

      MOAMu = this->levelShiftParam_ * MOAVir[mu];
      for(auto nu = 0; nu < this->nBasis_; nu++)
        FockA[nu + mu*this->nBasis_] += MOAMu * MOAVir[nu];

      if(this->nTCS_ == 1 && !this->isClosedShell){ 
        MOBMu = this->levelShiftParam_ * MOBVir[mu];
        for(auto nu = 0; nu < this->nBasis_; nu++)
          FockB[nu + mu*this->nBasis_] += MOBMu * MOBVir[nu];
      }
    
    }
  }
}
// JJRADLER
