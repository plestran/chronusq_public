/* 
 *  The Chronus Quantum (ChronusQ) software package is high-performace 
 *  computational chemistry software with a strong emphasis on explicitly 
 *  time-dependent and post-SCF quantum mechanical methods.
 *  
 *  Copyright (C) 2014-2015 Li Research Group (University of Washington)
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
  this->lenX_      = this->nBasis_ * this->nBasis_ * this->nTCS_ * this->nTCS_;
  this->lenXp_     = this->nBasis_ * this->nBasis_ * this->nTCS_ * this->nTCS_;
  this->lenF_      = this->nBasis_ * this->nBasis_ * this->nTCS_ * this->nTCS_;
  this->lenP_      = this->nBasis_ * this->nBasis_ * this->nTCS_ * this->nTCS_;
  this->lenCoeff_  = 7;
  this->lenB_      = this->lenCoeff_   * this->lenCoeff_;
  this->LWORK_     = 4 * this->nBasis_ * this->nTCS_;
  this->LRWORK_    = 3 * this->nBasis_ * this->nTCS_ - 2;
  this->lenLambda_ = this->nBasis_ * this->nBasis_ * this->nTCS_ * this->nTCS_;
  this->lenDelF_   = this->nBasis_ * this->nBasis_ * this->nTCS_ * this->nTCS_;
  this->lenOccNum_ = this->nBasis_ * this->nTCS_   * this->nTCS_;

  this->lenScr_     = 0;
  this->lenRealScr_ = 0;

  this->lenScr_ += this->lenX_;     // Storage for S^(-0.5)
  this->lenScr_ += this->lenF_;     // Storage for Alpha (Total) Fock
  this->lenScr_ += this->lenP_;     // Storage for Alpha (Total) Density
  this->lenScr_ += this->lenCoeff_; // Storage for CDIIS Coefficients
  this->lenScr_ += this->lenB_;     // Storage for CDIIS Metric
  this->lenScr_ += 2*(this->lenCoeff_ - 1) * this->lenF_; // CDIIS Commutator (A) array
  if(!this->isClosedShell && this->Ref_ != TCS) {
    this->lenScr_ += this->lenF_;     // Storage for Beta Fock
    this->lenScr_ += this->lenP_;     // Storage for Beta Density
    this->lenScr_ += 2*(this->lenCoeff_ - 1) * this->lenF_; // CDIIS Commutator (B) array
  } 
  if(this->Ref_ == CUHF) {
    this->lenRealScr_ += this->lenOccNum_; // Storage for Occupation Numbers (NOs)
    this->lenScr_     += this->lenXp_; // Storage for X^(0.5)
    this->lenScr_     += this->lenLambda_; // Storage for Lambda
    this->lenScr_     += this->lenDelF_;   // Stroage for DelF
    this->lenScr_     += this->lenP_;      // Storage for NOs
  }

  this->lenScr_ += this->LWORK_; // LAPACK Scratch space
  this->complexMem();

}; // initMemLen

template <typename T>
void SingleSlater<T>::initSCFPtr(){
  this->REAL_SCF_SCR   = NULL;
  this->occNumMem_     = NULL;
  this->RWORK_         = NULL;
  this->SCpyMem_       = NULL;
  this->SEVlMem_       = NULL;
  this->SEVcMem_       = NULL;
  this->LowdinWORK_    = NULL;

  this->SCF_SCR        = NULL;
  this->XMem_          = NULL;
  this->FpAlphaMem_    = NULL;
  this->FpBetaMem_     = NULL;
  this->POldAlphaMem_  = NULL;
  this->POldBetaMem_   = NULL;
  this->ErrorAlphaMem_ = NULL;
  this->ErrorBetaMem_  = NULL;
  this->FADIIS_        = NULL;
  this->FBDIIS_        = NULL;
  this->WORK_          = NULL;
  this->XpMem_         = NULL;
  this->lambdaMem_     = NULL;
  this->delFMem_       = NULL;
  this->PNOMem_        = NULL;
}; // initSCFPtr

template <typename T>
void SingleSlater<T>::initSCFMem(){
  this->initSCFPtr();
  this->initMemLen();

/*  This kind of allocation should be handled by a memory manager
  T* LAST_FOR_SECTION;
  int LEN_LAST_FOR_SECTION;

  this->SCF_SCR      = new T[this->lenScr_];
  this->REAL_SCF_SCR = new double[this->lenRealScr_];
  std::memset(this->SCF_SCR,     0.0,this->lenScr_    *sizeof(T)     );
  std::memset(this->REAL_SCF_SCR,0.0,this->lenRealScr_*sizeof(double));

  this->XMem_          = this->SCF_SCR;
  this->FpAlphaMem_    = this->XMem_          + this->lenX_;
  this->POldAlphaMem_  = this->FpAlphaMem_    + this->lenF_;
  this->ErrorAlphaMem_ = this->POldAlphaMem_  + this->lenP_;
  this->FADIIS_        = this->ErrorAlphaMem_ + this->lenF_*(this->lenCoeff_ -1);
  LAST_FOR_SECTION     = this->FADIIS_;
  LEN_LAST_FOR_SECTION = this->lenF_*(this->lenCoeff_ -1);
  if(!this->isClosedShell && this->Ref_ != TCS){
    this->FpBetaMem_     = LAST_FOR_SECTION + LEN_LAST_FOR_SECTION;
    this->POldBetaMem_   = this->FpBetaMem_    + this->lenF_;
    this->ErrorBetaMem_  = this->POldBetaMem_  + this->lenP_;
    this->FBDIIS_        = this->ErrorBetaMem_ + this->lenF_*(this->lenCoeff_ -1);
    LAST_FOR_SECTION     = this->FBDIIS_;
    LEN_LAST_FOR_SECTION = this->lenF_*(this->lenCoeff_ -1);
  }
  if(this->Ref_ == CUHF) {
    this->XpMem_     = LAST_FOR_SECTION + LEN_LAST_FOR_SECTION;
    this->delFMem_   = this->XpMem_     + this->lenX_;
    this->lambdaMem_ = this->delFMem_   + this->lenDelF_;
    this->PNOMem_    = this->lambdaMem_ + this->lenLambda_;
  //this->occNumMem_ = this->PNOMem_    + this->lenP_;
  //LAST_FOR_SECTION = this->occNumMem_;
  //LEN_LAST_FOR_SECTION = this->lenOccNum_;
    this->occNumMem_ = this->REAL_SCF_SCR;
    this->RWORK_     = this->occNumMem_ + this->lenOccNum_;
    LAST_FOR_SECTION = this->PNOMem_;
    LEN_LAST_FOR_SECTION = this->lenP_; 
  } else {
    this->RWORK_     = this->REAL_SCF_SCR;
  }
  
  this->WORK_  = LAST_FOR_SECTION + LEN_LAST_FOR_SECTION;
*/
  this->allocLowdin();
  this->allocAlphaScr();
  if(!this->isClosedShell && this->Ref_ != TCS) this->allocBetaScr();
  if(this->Ref_ == CUHF) this->allocCUHFScr();
  this->allocLAPACKScr();

  
}; //initSCFMem

template<typename T>
void SingleSlater<T>::allocLowdin(){
  auto NTCSxNBASIS = this->nBasis_*this->nTCS_;
  this->XMem_    = new T[this->lenX_];
  this->SCpyMem_ = new double[this->lenX_];
  this->SEVcMem_ = new double[this->lenX_];
  this->SEVlMem_ = new double[NTCSxNBASIS];
  if(typeid(T).hash_code() == typeid(dcomplex).hash_code()){
//  this->fileio_->out << "Allocating Extra LAPACK WORK in SCF" << endl;
    this->LowdinWORK_ = new double[this->LWORK_];
  }
}

template<typename T>
void SingleSlater<T>::allocAlphaScr(){
  this->FpAlphaMem_    = new T[this->lenF_];
  this->POldAlphaMem_  = new T[this->lenP_];
  this->ErrorAlphaMem_ = new T[this->lenF_*(this->lenCoeff_ -1)];
  this->FADIIS_        = new T[this->lenF_*(this->lenCoeff_ -1)];
}

template<typename T>
void SingleSlater<T>::allocBetaScr(){
  this->FpBetaMem_    = new T[this->lenF_];
  this->POldBetaMem_  = new T[this->lenP_];
  this->ErrorBetaMem_ = new T[this->lenF_*(this->lenCoeff_ -1)];
  this->FBDIIS_       = new T[this->lenF_*(this->lenCoeff_ -1)];
}

template<typename T>
void SingleSlater<T>::allocCUHFScr(){
  this->XpMem_     = new T[this->lenX_];
  this->delFMem_   = new T[this->lenDelF_];
  this->lambdaMem_ = new T[this->lenLambda_];
  this->PNOMem_    = new T[this->lenP_];
  this->occNumMem_ = new double[this->lenOccNum_];
}

template<typename T>
void SingleSlater<T>::allocLAPACKScr(){
  this->WORK_  = new T[this->LWORK_];
  if(typeid(T).hash_code() == typeid(dcomplex).hash_code()){
//  this->fileio_->out << "Allocating RWORK in SCF" << endl;
    this->RWORK_ = new double[this->LRWORK_];
  }
}

template<typename T>
void SingleSlater<T>::cleanupLowdin(){
  delete [] this->XMem_;
  delete [] this->SCpyMem_;
  delete [] this->SEVcMem_;
  delete [] this->SEVlMem_;
  delete [] this->LowdinWORK_;
}

template<typename T>
void SingleSlater<T>::cleanupSCFMem(){
  this->cleanupLowdin();
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
  delete [] this->XpMem_     ;
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
  if(this->printLevel_ > 0 && getRank() == 0) {
    this->printSCFHeader(this->fileio_->out);

    this->fileio_->out << std::setw(16) << "SCF Iteration";
    this->fileio_->out << std::setw(18) << "Energy (Eh)";
    this->fileio_->out << std::setw(18) << "\u0394E (Eh)";
    if(this->Ref_ == TCS)
      this->fileio_->out << std::setw(18) << "|\u0394P|";
    else {
      this->fileio_->out << std::setw(18) << "|\u0394P(\u03B1)|";
      if(!this->isClosedShell)
        this->fileio_->out << std::setw(18) << "|\u0394P(\u03B2)|";
    }
    this->fileio_->out << endl;
    this->fileio_->out << std::setw(16) << "-------------";
    this->fileio_->out << std::setw(18) << "-----------";
    this->fileio_->out << std::setw(18) << "-------";
    if(this->Ref_ == TCS)
      this->fileio_->out << std::setw(18) << "----";
    else {
      this->fileio_->out << std::setw(18) << "-------";
      if(!this->isClosedShell)
        this->fileio_->out << std::setw(18) << "-------";
    }
    this->fileio_->out << endl;
  }
  if(getRank() == 0) {
    this->initSCFMem();
    this->formX();
  }
  for (iter = 0; iter < this->maxSCFIter_; iter++){
/*
    this->fileio_->out << endl << endl << bannerTop <<endl;  
    this->fileio_->out << "SCF iteration:"<< iter+1 <<endl;  
    this->fileio_->out << bannerEnd <<endl;  
*/

    if(getRank() == 0) {
      if(this->Ref_ == CUHF) this->formNO();
      this->diagFock();
      if(iter == 0 && this->guess_ != READ) this->mixOrbitalsSCF();
    }
    this->formDensity();
    this->formFock();

    if(getRank() == 0) {
      if(this->Ref_ != CUHF && this->doDIIS){ // DIIS NYI for CUHF
        this->GenDComm(iter);
        this->CpyFock(iter);   
        if(iter % (this->lenCoeff_-1) == (this->lenCoeff_-2) && iter != 0) 
          this->CDIIS();
      }
    }

    this->evalConver(iter);
    this->nSCFIter++;
    if(this->isConverged) break;

  }; // SCF Loop
/*
  delete [] this->SCF_SCR;
  delete [] this->REAL_SCF_SCR;
*/
  if(getRank() == 0) this->cleanupSCFMem();
  printf("Hello 666%d:%d",getRank(),getSize());

  if(!this->isConverged && getRank() == 0)
    CErr("SCF Failed to converge within maximum number of iterations",this->fileio_->out);

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

