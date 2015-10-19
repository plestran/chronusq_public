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
void RealTime<T>::iniRealTime(FileIO *fileio, Controls *controls, 
       AOIntegrals *aointegrals, SingleSlater<T> *groundState) {
  this->communicate(*fileio,*controls,*aointegrals,*groundState);
  this->initMeta();


// Keep around for C++ interface (for now)
  this->maxSteps_	= this->controls_->rtMaxSteps;
  this->stepSize_	= this->controls_->rtTimeStep;
  this->typeOrtho_	= this->controls_->rtTypeOrtho; 
  this->initDensity_	= this->controls_->rtInitDensity;
  this->swapMOA_	= this->controls_->rtSwapMOA;
  this->swapMOB_	= this->controls_->rtSwapMOB;
  this->methFormU_	= this->controls_->rtMethFormU;
  this->IEnvlp_         = this->controls_->rtEnvelope_;
  this->Ex_             = this->controls_->rtField_[0];
  this->Ey_             = this->controls_->rtField_[1];
  this->Ez_             = this->controls_->rtField_[2];
  this->TOn_            = this->controls_->rtTOn_;
  this->TOff_           = this->controls_->rtTOff_;
  this->Freq_          = this->controls_->rtFreq_;
  this->Phase_          = this->controls_->rtPhase_;
  this->Sigma_          = this->controls_->rtSigma_;
  this->printLevel_     = this->controls_->printLevel;

  this->alloc();
};

template<typename T>
void RealTime<T>::alloc(){
  this->checkMeta();

  this->ssPropagator_	= std::unique_ptr<SingleSlater<dcomplex>>(
                            new SingleSlater<dcomplex>(
                              this->groundState_));
  this->initMem();
}

template<typename T>
void RealTime<T>::initRTPtr(){
  this->SCR         = NULL;
  this->oTrans1Mem_ = NULL;
  this->oTrans2Mem_ = NULL;
  this->POAMem_     = NULL;
  this->POBMem_     = NULL;
  this->POAsavMem_  = NULL;
  this->POBsavMem_  = NULL;
  this->FOAMem_     = NULL;
  this->FOBMem_     = NULL;
  this->initMOAMem_ = NULL;
  this->initMOBMem_ = NULL;
  this->uTransAMem_ = NULL;
  this->uTransBMem_ = NULL;
  this->scratchMem_ = NULL;

  this->REAL_LAPACK_SCR  = NULL;
  this->CMPLX_LAPACK_SCR = NULL;
}

template<typename T>
void RealTime<T>::initMemLen(){
  auto NTCSxNBASIS = this->nTCS_*this->nBasis_;

  this->lenOTrans1_ = NTCSxNBASIS * NTCSxNBASIS;
  this->lenOTrans2_ = NTCSxNBASIS * NTCSxNBASIS;
  this->lenPOA_     = NTCSxNBASIS * NTCSxNBASIS;
  this->lenPOB_     = NTCSxNBASIS * NTCSxNBASIS;
  this->lenPOAsav_  = NTCSxNBASIS * NTCSxNBASIS;
  this->lenPOBsav_  = NTCSxNBASIS * NTCSxNBASIS;
  this->lenFOA_     = NTCSxNBASIS * NTCSxNBASIS;
  this->lenFOB_     = NTCSxNBASIS * NTCSxNBASIS;
  this->lenInitMOA_ = NTCSxNBASIS * NTCSxNBASIS;
  this->lenInitMOB_ = NTCSxNBASIS * NTCSxNBASIS;
  this->lenUTransA_ = NTCSxNBASIS * NTCSxNBASIS;
  this->lenUTransB_ = NTCSxNBASIS * NTCSxNBASIS;
  this->lenScratch_ = NTCSxNBASIS * NTCSxNBASIS;

  this->lenScr_ = 0;
  
  this->lenScr_ += this->lenOTrans1_; 
  this->lenScr_ += this->lenOTrans2_; 
  this->lenScr_ += this->lenPOA_    ; 
  this->lenScr_ += this->lenPOAsav_ ; 
  this->lenScr_ += this->lenFOA_    ; 
  this->lenScr_ += this->lenInitMOA_; 
  this->lenScr_ += this->lenUTransA_; 
  this->lenScr_ += this->lenScratch_; 
  if(!this->isClosedShell_ && this->Ref_ != SingleSlater<T>::TCS){
    this->lenScr_ += this->lenPOB_    ; 
    this->lenScr_ += this->lenPOBsav_ ; 
    this->lenScr_ += this->lenFOB_    ; 
    this->lenScr_ += this->lenInitMOB_; 
    this->lenScr_ += this->lenUTransB_; 
  }

  this->lWORK               =  NTCSxNBASIS * NTCSxNBASIS;
  this->lenREAL_LAPACK_SCR  =  0;
  this->lenREAL_LAPACK_SCR  += NTCSxNBASIS * NTCSxNBASIS;   // DSYEV::A
  this->lenREAL_LAPACK_SCR  += NTCSxNBASIS;                 // DSYEV/ZHEEV::W
  this->lenREAL_LAPACK_SCR  += this->lWORK;                 // DSYEV::WORK
  this->lenREAL_LAPACK_SCR  += std::max(1,3*NTCSxNBASIS-2); // ZHEEV::RWORK
  this->lenCMPLX_LAPACK_SCR =  0;
  // This can be handled my the scratch space
//this->lenCMPLX_LAPACK_SCR += NTCSxNBASIS * NTCSxNBASIS;   // ZHEEV::A
  this->lenCMPLX_LAPACK_SCR += this->lWORK;                 // ZHEEV::WORK

  this->lenScr_ += this->lenCMPLX_LAPACK_SCR;
}

template<typename T>
void RealTime<T>::initMem(){
  try{ this->SCR = new dcomplex[this->lenScr_]; } 
  catch (...) { CErr(std::current_exception(),"RT Allocation"); }
  dcomplex *LAST_OF_SECTION ;
  int  LEN_LAST_OF_SECTION  ;

  this->oTrans1Mem_ = this->SCR;
  this->oTrans2Mem_ = this->oTrans1Mem_ + this->lenOTrans1_;
  this->POAMem_     = this->oTrans2Mem_ + this->lenOTrans2_;
  this->POAsavMem_  = this->POAMem_     + this->lenPOA_;
  this->FOAMem_     = this->POAsavMem_  + this->lenPOAsav_;
  this->initMOAMem_ = this->FOAMem_     + this->lenFOA_;
  this->uTransAMem_ = this->initMOAMem_ + this->lenInitMOA_;
  LAST_OF_SECTION      = this->uTransAMem_;
  LEN_LAST_OF_SECTION  = this->lenUTransA_;

  if(!this->isClosedShell_ && this->Ref_ != SingleSlater<T>::TCS){
    this->POBMem_     = LAST_OF_SECTION   + LEN_LAST_OF_SECTION;
    this->POBsavMem_  = this->POBMem_     + this->lenPOB_;
    this->FOBMem_     = this->POBsavMem_  + this->lenPOBsav_;
    this->initMOBMem_ = this->FOBMem_     + this->lenFOB_;
    this->uTransBMem_ = this->initMOBMem_ + this->lenInitMOB_;
    LAST_OF_SECTION      = this->uTransBMem_;
    LEN_LAST_OF_SECTION  = this->lenUTransB_;
  }

  this->scratchMem_      = LAST_OF_SECTION       + LEN_LAST_OF_SECTION;
  this->CMPLX_LAPACK_SCR = this->scratchMem_     + this->lenScratch_;

  this->REAL_LAPACK_SCR = new double[this->lenREAL_LAPACK_SCR];
}

