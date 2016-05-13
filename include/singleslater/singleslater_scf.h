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
void SingleSlater<T>::initMemLen(){
  this->lenF_      = this->nBasis_ * this->nBasis_ * this->nTCS_ * this->nTCS_;
  this->lenP_      = this->nBasis_ * this->nBasis_ * this->nTCS_ * this->nTCS_;
  this->lenCoeff_  = this->nDIISExtrap_;
  this->lenB_      = this->lenCoeff_   * this->lenCoeff_;
  this->LWORK_     = 5 * std::max(this->nBasis_ * this->nTCS_,
    this->nDIISExtrap_);
  this->LRWORK_    = 3 * 
    std::max(this->nBasis_ * this->nTCS_,this->nDIISExtrap_) - 2;
  this->lenLambda_ = this->nBasis_ * this->nBasis_ * this->nTCS_ * this->nTCS_;
  this->lenDelF_   = this->nBasis_ * this->nBasis_ * this->nTCS_ * this->nTCS_;
  this->lenOccNum_ = this->nBasis_ * this->nTCS_   * this->nTCS_;
}; // initMemLen

template <typename T>
void SingleSlater<T>::initSCFPtr(){
  this->occNumMem_     = NULL;
  this->RWORK_         = NULL;

  this->FpAlphaMem_    = NULL;
  this->FpBetaMem_     = NULL;
  this->POldAlphaMem_  = NULL;
  this->POldBetaMem_   = NULL;
  this->ErrorAlphaMem_ = NULL;
  this->ErrorBetaMem_  = NULL;
  this->FADIIS_        = NULL;
  this->FBDIIS_        = NULL;
  this->WORK_          = NULL;
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

template <typename T>
void SingleSlater<T>::initSCFMem(){
  this->initSCFPtr();
  this->initMemLen();

  this->allocAlphaScr();
  if(!this->isClosedShell && this->Ref_ != TCS) this->allocBetaScr();
  if(this->Ref_ == CUHF) this->allocCUHFScr();
  this->allocLAPACKScr();

  
}; //initSCFMem


template<typename T>
void SingleSlater<T>::allocAlphaScr(){
  this->FpAlphaMem_    = new T[this->lenF_];
  this->POldAlphaMem_  = new T[this->lenP_];
  this->ErrorAlphaMem_ = new T[this->lenF_*(this->nDIISExtrap_ -1)];
  this->FADIIS_        = new T[this->lenF_*(this->nDIISExtrap_ -1)];
}

template<typename T>
void SingleSlater<T>::allocBetaScr(){
  this->FpBetaMem_    = new T[this->lenF_];
  this->POldBetaMem_  = new T[this->lenP_];
  this->ErrorBetaMem_ = new T[this->lenF_*(this->nDIISExtrap_ -1)];
  this->FBDIIS_       = new T[this->lenF_*(this->nDIISExtrap_ -1)];
}

template<typename T>
void SingleSlater<T>::allocCUHFScr(){
  this->delFMem_   = new T[this->lenDelF_];
  this->lambdaMem_ = new T[this->lenLambda_];
  this->PNOMem_    = new T[this->lenP_];
  this->occNumMem_ = new double[this->lenOccNum_];
}

template<typename T>
void SingleSlater<T>::allocLAPACKScr(){
  this->WORK_  = new T[this->LWORK_];
  if(typeid(T).hash_code() == typeid(dcomplex).hash_code()){
    this->RWORK_ = new double[this->LRWORK_];
  }
}


template<typename T>
void SingleSlater<T>::cleanupSCFMem(){
  //this->cleanupLowdin();
  this->cleanupAlphaScr();
  if(!this->isClosedShell && this->Ref_ != TCS) this->cleanupBetaScr();
  if(this->Ref_ == CUHF) this->cleanupCUHFScr();
  this->cleanupLAPACKScr();
}

template<typename T>
void SingleSlater<T>::cleanupAlphaScr(){
  delete [] this->FpAlphaMem_    ;
  delete [] this->POldAlphaMem_  ;
  delete [] this->ErrorAlphaMem_ ;
  delete [] this->FADIIS_        ;
}

template<typename T>
void SingleSlater<T>::cleanupBetaScr(){
  delete [] this->FpBetaMem_    ;
  delete [] this->POldBetaMem_  ;
  delete [] this->ErrorBetaMem_ ;
  delete [] this->FBDIIS_       ;
}

template<typename T>
void SingleSlater<T>::cleanupCUHFScr(){
  delete [] this->delFMem_   ;
  delete [] this->lambdaMem_ ;
  delete [] this->PNOMem_    ;
  delete [] this->occNumMem_ ;
}

template<typename T>
void SingleSlater<T>::cleanupLAPACKScr(){
  delete [] this->WORK_  ;
  delete [] this->RWORK_ ;
}

template<typename T>
void SingleSlater<T>::SCF(){

  if(!this->aointegrals_->haveAOOneE && getRank() == 0) 
    this->aointegrals_->computeAOOneE();

  int iter; 
  if(getRank() == 0) {
    if(this->printLevel_ > 0)
      this->printSCFHeader(this->fileio_->out);
    this->initSCFMem();
  }
  for (iter = 0; iter < this->maxSCFIter_; iter++){

    if(getRank() == 0) {
      if(this->Ref_ == CUHF) this->formNO();
      this->diagFock();
      if(iter == 0 && this->guess_ != READ) this->mixOrbitalsSCF();
    }
    this->formDensity();
    this->formFock();
    if(PyErr_CheckSignals() == -1)
      CErr("Keyboard Interrupt in SCF!",this->fileio_->out);
   

    if(getRank() == 0) {
      if(iter == this->iDIISStart_ && this->printLevel_ > 0)
        this->fileio_->out << 
	  std::setw(2) << " " <<
	  std::setw(4) << " " <<
	  "*** Starting DIIS ***" << endl;
      // DIIS NYI for CUHF
      if(this->Ref_ != CUHF && this->doDIIS && iter >= this->iDIISStart_){ 
        this->GenDComm(iter - this->iDIISStart_);
        this->CpyFock(iter - this->iDIISStart_);   
        if((iter - this->iDIISStart_) % (this->nDIISExtrap_-1) == 
	   (this->nDIISExtrap_-2) && iter != 0) 
          this->CDIIS();
      }
    }

    this->evalConver(iter);
    this->nSCFIter++;
    if(this->isConverged) break;

  }; // SCF Loop
//prettyPrint(cout,(*this->aointegrals_->overlap_),"MOS");
/*
  delete [] this->SCF_SCR;
  delete [] this->REAL_SCF_SCR;
*/
  if(getRank() == 0){
    this->cleanupSCFMem();
    this->fixPhase();
  }

  if(!this->isConverged && getRank() == 0)
    CErr("SCF Failed to converge within maximum number of iterations",
        this->fileio_->out);

  if(this->printLevel_ > 0 && getRank() == 0){
    if(this->isConverged){
      this->fileio_->out << endl << "SCF Completed: E(" << this->SCFTypeShort_ 
                         << ") = ";
      this->fileio_->out << std::fixed << std::setprecision(10) 
                         << this->totalEnergy << "  Eh after  " << iter + 1 
                         << "  SCF Iterations" << endl;
    }
    this->fileio_->out << bannerEnd <<endl;
  }
}


template<typename T>
void SingleSlater<T>::initSCFMem2(){
  this->initSCFPtr();


  auto NTCSxNBASIS = this->nTCS_ * this->nBasis_;
  auto NSQ = NTCSxNBASIS  * NTCSxNBASIS;

  this->LWORK_     = 5 * std::max(NTCSxNBASIS,this->nDIISExtrap_);
  this->LRWORK_    = 3 * std::max(NTCSxNBASIS,this->nDIISExtrap_);
  this->lenCoeff_  = this->nDIISExtrap_;

  /*
//  this->FpAlphaMem_    = new T[NSQ];
  this->POldAlphaMem_  = new T[NSQ];
  this->ErrorAlphaMem_ = new T[NSQ*(this->nDIISExtrap_ -1)];
  this->FADIIS_        = new T[NSQ*(this->nDIISExtrap_ -1)];
  if(this->nTCS_ == 1 && !this->isClosedShell) {
//    this->FpBetaMem_    = new T[NSQ*NTCSxNBASIS];
    this->POldBetaMem_  = new T[NSQ*NTCSxNBASIS];
    this->ErrorBetaMem_ = new T[NSQ*NTCSxNBASIS*(this->nDIISExtrap_ -1)];
    this->FBDIIS_       = new T[NSQ*NTCSxNBASIS*(this->nDIISExtrap_ -1)];
  }


  if(this->Ref_ == CUHF) {
    this->delFMem_   = new T[NSQ];
    this->lambdaMem_ = new T[NSQ];
    this->PNOMem_    = new T[NSQ];
    this->occNumMem_ = new double[NTCSxNBASIS];
  }

  this->WORK_  = new T[this->LWORK_];
  if(typeid(T).hash_code() == typeid(dcomplex).hash_code())
    this->RWORK_ = new double[this->LRWORK_];
  */
  
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

  this->WORK_  = this->memManager_->template malloc<T>(this->LWORK_);
  if(typeid(T).hash_code() == typeid(dcomplex).hash_code())
    this->RWORK_ = this->memManager_->template malloc<double>(this->LRWORK_);
  
}; //initSCFMem

template<typename T>
void SingleSlater<T>::SCF2(){
  this->aointegrals_->computeAOOneE();
  if(this->printLevel_ > 0)
    this->printSCFHeader(this->fileio_->out);
  this->initSCFMem2();

  size_t iter;
  for(iter = 0; iter < this->maxSCFIter_; iter++){
    if(this->Ref_ == CUHF) {
      this->formNO();
      this->fockCUHF();
    }

    this->orthoFock();
    this->diagFock2();
    if(iter == 0 && this->guess_ != READ) this->mixOrbitalsSCF();

    this->copyDen();
    this->formDensity();
    this->orthoDen();
    this->formFock();

    if(PyErr_CheckSignals() == -1)
      CErr("Keyboard Interrupt in SCF!",this->fileio_->out);
    if(iter == this->iDIISStart_ && this->printLevel_ > 0)
      this->fileio_->out << 
        std::setw(2) << " " <<
        std::setw(4) << " " <<
        "*** Starting DIIS ***" << endl;
    // DIIS NYI for CUHF
    if(this->Ref_ != CUHF && this->doDIIS && iter >= this->iDIISStart_){ 
      this->genDComm2(iter - this->iDIISStart_);
      this->CpyFock(iter - this->iDIISStart_);   
      if((iter - this->iDIISStart_) % (this->nDIISExtrap_-1) == 
         (this->nDIISExtrap_-2) && iter != 0) 
        this->CDIIS();
    }

    this->evalConver(iter);
    this->nSCFIter++;

    if(this->isConverged) break;
  };

  this->cleanupSCFMem2();
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
  /*
//  delete[] this->FpAlphaMem_    
  delete[] this->POldAlphaMem_;  
  delete[] this->ErrorAlphaMem_; 
  delete[] this->FADIIS_;        
  if(this->nTCS_ == 2 && !this->isClosedShell) { 
//    delete[] this->FpBetaMem_;    
    delete[] this->POldBetaMem_;  
    delete[] this->ErrorBetaMem_; 
    delete[] this->FBDIIS_;       
  }

  if(this->Ref_ == CUHF){ 
    delete[] this->delFMem_;   
    delete[] this->lambdaMem_; 
    delete[] this->PNOMem_;    
    delete[] this->occNumMem_; 
  }

  delete[] this->WORK_;  
  if(typeid(T).hash_code() == typeid(dcomplex).hash_code()) 
    delete[] this->RWORK_; 
    */
                                  
  auto NTCSxNBASIS = this->nTCS_ * this->nBasis_;
  auto NSQ = NTCSxNBASIS  * NTCSxNBASIS;

  this->memManager_->free(this->POldAlphaMem_,NSQ);  
  this->memManager_->free(this->ErrorAlphaMem_,(this->nDIISExtrap_-1)*NSQ); 
  this->memManager_->free(this->FADIIS_,(this->nDIISExtrap_-1)*NSQ);        
  if(this->nTCS_ == 2 && !this->isClosedShell) { 
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

  this->memManager_->free(this->WORK_,this->LWORK_);  
  if(typeid(T).hash_code() == typeid(dcomplex).hash_code()) 
    this->memManager_->free(this->RWORK_,this->LRWORK_); 
  
  
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
