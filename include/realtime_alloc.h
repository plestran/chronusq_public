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

  this->fileio_->out << endl << "Real-time TDHF: "<< endl;
  this->fileio_->out << std::right << std::setw(20) << "Number of steps = "
                     << std::setw(15) << this->maxSteps_ << std::setw(5)
                     << endl;
  this->fileio_->out << std::right << std::setw(20) << "Step size = "
                     <<std::setw(15) << this->stepSize_ << std::setw(5) 
                     << " a.u. " << endl;
  this->alloc();
};

template<typename T>
void RealTime<T>::alloc(){
  this->checkMeta();

  this->ssPropagator_	= std::unique_ptr<SingleSlater<dcomplex>>(
                            new SingleSlater<dcomplex>(
                              this->groundState_));

  auto NTCSxNBASIS = this->nTCS_*this->nBasis_;
  try {
    this->oTrans1_ = 
      std::unique_ptr<ComplexMatrix>(
        new ComplexMatrix(NTCSxNBASIS,NTCSxNBASIS));
  } catch(...) {
    CErr(std::current_exception(),"Orthog. Trans. Matrix 1");
  }
  try {
    this->oTrans2_ = 
      std::unique_ptr<ComplexMatrix>(
        new ComplexMatrix(NTCSxNBASIS,NTCSxNBASIS));
  } catch(...) {
    CErr(std::current_exception(),"Orthog. Trans. Matrix 2 (Inverse)");
  }
  try {
    this->POA_  	 = 
      std::unique_ptr<ComplexMatrix>(
        new ComplexMatrix(NTCSxNBASIS,NTCSxNBASIS));
  } catch(...) {
    CErr(std::current_exception(),"Density Alpha in Orthonormal Basis");
  }
  try {
    this->POAsav_  = 
      std::unique_ptr<ComplexMatrix>(
        new ComplexMatrix(NTCSxNBASIS,NTCSxNBASIS));
  } catch(...) {
    CErr(std::current_exception(),"copy of Density Alpha in Orthonormal Basis");
  }
  try {
    this->FOA_ 	 = 
      std::unique_ptr<ComplexMatrix>(
        new ComplexMatrix(NTCSxNBASIS,NTCSxNBASIS));
  } catch(...) {
    CErr(std::current_exception(),"Fock Matrix Alpha in Ortho. Basis");
  }
  try {
    this->initMOA_ = 
      std::unique_ptr<ComplexMatrix>(
        new ComplexMatrix(NTCSxNBASIS,NTCSxNBASIS));
  } catch(...) {
    CErr(std::current_exception(),"GS MO Alpha in orthonormal basis");
  }
  try {
    this->uTransA_ = 
      std::unique_ptr<ComplexMatrix>(
        new ComplexMatrix(NTCSxNBASIS,NTCSxNBASIS));
  } catch(...) {
    CErr(std::current_exception(),"Unitary Trans. Matrix Alpha [exp(-i*dt*F)] ");
  }
  try {
    this->scratch_ = 
      std::unique_ptr<ComplexMatrix>(
        new ComplexMatrix(NTCSxNBASIS,NTCSxNBASIS));
  } catch(...) {
    CErr(std::current_exception(),"N x N scratch for RT");
  }

  if(!this->isClosedShell_ && this->Ref_ != SingleSlater<T>::TCS){
    try {
      this->POB_  	 = 
        std::unique_ptr<ComplexMatrix>(
          new ComplexMatrix(NTCSxNBASIS,NTCSxNBASIS));
    } catch(...) {
      CErr(std::current_exception(),"Density Beta in Ortho. Basis");
    }
    try {
      this->POBsav_  = 
        std::unique_ptr<ComplexMatrix>(
          new ComplexMatrix(NTCSxNBASIS,NTCSxNBASIS));
    } catch(...) {
      CErr(std::current_exception(),"copy of Density Beta in Ortho. Basis");
    }
    try {
      this->FOB_ 	 = 
        std::unique_ptr<ComplexMatrix>(
          new ComplexMatrix(NTCSxNBASIS,NTCSxNBASIS));
    } catch(...) {
      CErr(std::current_exception(),"Fock Matrix Beta in Ortho. Basis");
    }
    try {
      this->initMOB_ = 
        std::unique_ptr<ComplexMatrix>(
          new ComplexMatrix(NTCSxNBASIS,NTCSxNBASIS));
    } catch(...) {
      CErr(std::current_exception(),"GS MO Beta in Ortho. Basis");
    }
    try {
      this->uTransB_ = 
        std::unique_ptr<ComplexMatrix>(
          new ComplexMatrix(NTCSxNBASIS,NTCSxNBASIS));
    } catch(...) {
      CErr(std::current_exception(),"Unitary Trans. Matrix Beta [exp(-i*dt*F)]");
    }
  }
}
