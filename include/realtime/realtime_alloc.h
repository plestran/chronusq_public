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
  this->Freq_           = this->controls_->rtFreq_;
  this->Phase_          = this->controls_->rtPhase_;
  this->Sigma_          = this->controls_->rtSigma_;
  this->printLevel_     = this->controls_->printLevel;

  this->alloc();
};

template<typename T>
void RealTime<T>::alloc(){
  this->checkMeta();

//cout << "HERE" << endl;
  this->ssPropagator_	= std::unique_ptr<SingleSlater<dcomplex>>(
                            new SingleSlater<dcomplex>(
                              this->groundState_));
  this->ssPropagator_->isNotPrimary();
//prettyPrintComplex(cout,*this->ssPropagator_->onePDMA(),"PA In RT Alloc1");
  if(getRank() == 0) this->initMem();
//prettyPrintComplex(cout,*this->ssPropagator_->onePDMA(),"PA In RT Alloc2");
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
  this->scratchMem2_ = NULL;

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
  this->lenScratch2_ = NTCSxNBASIS * NTCSxNBASIS;

  this->lenScr_ = 0;
  
  this->lenScr_ += this->lenOTrans1_; 
  this->lenScr_ += this->lenOTrans2_; 
  this->lenScr_ += this->lenPOA_    ; 
  this->lenScr_ += this->lenPOAsav_ ; 
  this->lenScr_ += this->lenFOA_    ; 
  this->lenScr_ += this->lenInitMOA_; 
  this->lenScr_ += this->lenUTransA_; 
  this->lenScr_ += this->lenScratch_; 
  this->lenScr_ += this->lenScratch2_; 
  if(!this->isClosedShell_ && this->Ref_ != SingleSlater<T>::TCS){
    this->lenScr_ += this->lenPOB_    ; 
    this->lenScr_ += this->lenPOBsav_ ; 
    this->lenScr_ += this->lenFOB_    ; 
    this->lenScr_ += this->lenInitMOB_; 
    this->lenScr_ += this->lenUTransB_; 
  }

  this->lWORK               =  NTCSxNBASIS*NTCSxNBASIS + 1;   // DSYEV/ZHEEV::LWORK 
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

  this->oTrans1Mem_ = 
    this->memManager_->template malloc<dcomplex>(this->lenOTrans1_);
  this->oTrans2Mem_ = 
    this->memManager_->template malloc<dcomplex>(this->lenOTrans2_);
  this->POAMem_     = 
    this->memManager_->template malloc<dcomplex>(this->lenPOA_);
  this->POAsavMem_  = 
    this->memManager_->template malloc<dcomplex>(this->lenPOAsav_);
  this->FOAMem_     = 
    this->memManager_->template malloc<dcomplex>(this->lenFOA_);
  this->initMOAMem_  = 
    this->memManager_->template malloc<dcomplex>(this->lenInitMOA_);
  this->uTransAMem_ = 
    this->memManager_->template malloc<dcomplex>(this->lenUTransA_);

  if(this->nTCS_ == 1 && !this->isClosedShell_){
    this->POBMem_     = 
      this->memManager_->template malloc<dcomplex>(this->lenPOB_);
    this->POBsavMem_  = 
      this->memManager_->template malloc<dcomplex>(this->lenPOBsav_);
    this->FOBMem_     = 
      this->memManager_->template malloc<dcomplex>(this->lenFOB_);
    this->initMOBMem_  = 
      this->memManager_->template malloc<dcomplex>(this->lenInitMOB_);
    this->uTransBMem_ = 
      this->memManager_->template malloc<dcomplex>(this->lenUTransB_);
  }

  this->scratchMem_  = 
    this->memManager_->template malloc<dcomplex>(this->lenScratch_);
  this->scratchMem2_ = 
    this->memManager_->template malloc<dcomplex>(this->lenScratch2_);

  this->CMPLX_LAPACK_SCR = 
    this->memManager_->template malloc<dcomplex>(this->lenCMPLX_LAPACK_SCR);
  this->REAL_LAPACK_SCR =
    this->memManager_->template malloc<double>(this->lenREAL_LAPACK_SCR);

  std::fill_n(this->oTrans1Mem_,this->lenOTrans1_, dcomplex(0.0,0.0));
  std::fill_n(this->oTrans2Mem_,this->lenOTrans2_, dcomplex(0.0,0.0));
  std::fill_n(this->POAMem_    ,this->lenPOA_,     dcomplex(0.0,0.0));
  std::fill_n(this->POAsavMem_ ,this->lenPOAsav_,  dcomplex(0.0,0.0));
  std::fill_n(this->FOAMem_    ,this->lenFOA_,     dcomplex(0.0,0.0));
  std::fill_n(this->initMOAMem_,this->lenInitMOA_, dcomplex(0.0,0.0));
  std::fill_n(this->uTransAMem_,this->lenUTransA_, dcomplex(0.0,0.0));
  if(this->nTCS_ == 1 && !this->isClosedShell_){
    std::fill_n(this->POBMem_    , this->lenPOB_,     dcomplex(0.0,0.0));
    std::fill_n(this->POBsavMem_ , this->lenPOBsav_,  dcomplex(0.0,0.0));
    std::fill_n(this->FOBMem_    , this->lenFOB_,     dcomplex(0.0,0.0));
    std::fill_n(this->initMOBMem_, this->lenInitMOB_, dcomplex(0.0,0.0));
    std::fill_n(this->uTransBMem_, this->lenUTransB_, dcomplex(0.0,0.0));
  }
  std::fill_n(this->scratchMem_ , this->lenScratch_,  dcomplex(0.0,0.0));
  std::fill_n(this->scratchMem2_, this->lenScratch2_, dcomplex(0.0,0.0));

  std::fill_n(this->CMPLX_LAPACK_SCR,this->lenCMPLX_LAPACK_SCR,
      dcomplex(0.0,0.0));
  std::fill_n(this->REAL_LAPACK_SCR ,this->lenREAL_LAPACK_SCR,0.0);
}

template<typename T>
void RealTime<T>::initCSV(){
  // Create/open CSVs for printing results
  csvs.push_back(
    new std::ofstream(this->fileio_->fileName() + "_RealTime_Dipole.csv")
  );
  this->csvFiles[this->csvs[0]] = 
    this->fileio_->fileName() + "_RealTime_Dipole.csv";

  csvs.push_back(
    new std::ofstream(this->fileio_->fileName() + "_RealTime_AppliedField.csv")
  );
  this->csvFiles[this->csvs[1]] = 
    this->fileio_->fileName() + "_RealTime_AppliedField.csv";

  csvs.push_back(
    new std::ofstream(this->fileio_->fileName() + "_RealTime_Mulliken.csv")
  );
  this->csvFiles[this->csvs[2]] = 
    this->fileio_->fileName() + "_RealTime_Mulliken.csv";

  csvs.push_back(
    new std::ofstream(this->fileio_->fileName() + "_RealTime_OrbOcc_Alpha.csv")
  );
  this->csvFiles[this->csvs[3]] = 
    this->fileio_->fileName() + "_RealTime_OrbOcc_Alpha.csv";

  if(!this->isClosedShell_ && this->Ref_ != SingleSlater<double>::TCS){
    csvs.push_back(
      new std::ofstream(this->fileio_->fileName() + "_RealTime_OrbOcc_Beta.csv")
    );
    this->csvFiles[this->csvs[4]] = 
      this->fileio_->fileName() + "_RealTime_OrbOcc_Beta.csv";
  }
};

template<typename T>
void RealTime<T>::cleanup() {
  this->memManager_->free(this->oTrans1Mem_ , this->lenOTrans1_);
  this->memManager_->free(this->oTrans2Mem_ , this->lenOTrans2_);
  this->memManager_->free(this->POAMem_     , this->lenPOA_);
  this->memManager_->free(this->POAsavMem_  , this->lenPOAsav_);
  this->memManager_->free(this->FOAMem_     , this->lenFOA_);
  this->memManager_->free(this->initMOAMem_  , this->lenInitMOA_);
  this->memManager_->free(this->uTransAMem_ , this->lenUTransA_);
  if(this->nTCS_ == 1 && !this->isClosedShell_){
    this->memManager_->free(this->POBMem_     , this->lenPOB_);
    this->memManager_->free(this->POBsavMem_  , this->lenPOBsav_);
    this->memManager_->free(this->FOBMem_     , this->lenFOB_);
    this->memManager_->free(this->initMOBMem_  , this->lenInitMOB_);
    this->memManager_->free(this->uTransBMem_ , this->lenUTransB_);
  }
  this->memManager_->free(this->scratchMem_  , this->lenScratch_);
  this->memManager_->free(this->scratchMem2_ , this->lenScratch2_);

  this->memManager_->free(this->CMPLX_LAPACK_SCR,this->lenCMPLX_LAPACK_SCR);
  this->memManager_->free(this->REAL_LAPACK_SCR,this->lenREAL_LAPACK_SCR);

  if(this->tarCSVs) this->tarCSVFiles();
}
