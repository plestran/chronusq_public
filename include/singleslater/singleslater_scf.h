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

#include <singleslater/singleslater_oldscf.h>
#include <singleslater/singleslater_levelshift.h>
#include <singleslater/singleslater_scfutils.h>


template <typename T>
void SingleSlater<T>::copyDen(){
/*
  auto NTCSxNBASIS = this->nTCS_*this->nBasis_;
  TMap POldAlpha(this->POldAlphaMem_,NTCSxNBASIS,NTCSxNBASIS);
  POldAlpha = (*this->onePDMA_);

  if(this->nTCS_ == 1 && !this->isClosedShell){
    TMap POldBeta(this->POldBetaMem_,NTCSxNBASIS,NTCSxNBASIS);
    POldBeta = (*this->onePDMB_);
  };
*/
  
  T* ScalarPtr;
  if(this->isClosedShell && this->nTCS_ == 1) 
    ScalarPtr = this->onePDMA_->data();
  else
    ScalarPtr = this->onePDMScalar_->data();

  DScalarOld_->write(ScalarPtr,H5PredType<T>());
  if(this->nTCS_ == 2 || !this->isClosedShell){
    DMzOld_->write(this->onePDMMz_->data(),H5PredType<T>());
  }
  if(this->nTCS_ == 2){
    DMyOld_->write(this->onePDMMy_->data(),H5PredType<T>());
    DMxOld_->write(this->onePDMMx_->data(),H5PredType<T>());
  }
};

template<typename T>
void SingleSlater<T>::SCF3(){
  // Compute the 1-Body Integrals if they haven't already been 
  // computed
  if(!this->aointegrals_->haveAOOneE) 
    this->aointegrals_->computeAOOneE();

  this->initSCFMem3();

  // Print the SCF Header
  if(this->printLevel_ > 0)
    this->printSCFHeader(this->fileio_->out);

  this->isConverged = false;
  std::size_t iter;
  this->fileio_->out << "    *** INITIAL GUESS ENERGY = " 
    << this->totalEnergy << " Eh ***" << endl;

  this->formFock();
  for(iter = 0; iter < this->maxSCFIter_; iter++){
    this->copyDen();
    this->orthoFock3();
    if(this->doDIIS){
      auto IDIISIter = iter - this->iDIISStart_;
      this->cpyFockDIIS(IDIISIter);
      this->cpyDenDIIS(IDIISIter);
      this->genDIISCom(IDIISIter); 
      if(IDIISIter > 0) 
        this->CDIIS4(std::min(IDIISIter+1,std::size_t(this->nDIISExtrap_)));
    }
    this->diagFock2();
    this->formDensity();
    this->cpyAOtoOrthoDen();
    this->unOrthoDen3();
    // Calls formFock
    SCFConvergence CONVER = this->evalConver3();

    if(this->printLevel_ > 0)
      this->printSCFIter(iter,CONVER.EDelta,CONVER.PARMS,
        CONVER.PBRMS);

    this->nSCFIter++;

    if(this->isConverged) break;
  }

  this->cleanupSCFMem3();
  this->backTransformMOs();
  this->fixPhase();
  if(!this->isConverged)
    CErr("SCF Failed to converge within MAXITER iterations",
        this->fileio_->out);
  if(this->printLevel_ > 0 ){
    this->fileio_->out 
      << endl << "SCF Completed: E(" 
      << this->SCFTypeShort_ << ") = ";
    this->fileio_->out 
      << std::fixed << std::setprecision(10) 
      << this->totalEnergy << "  Eh after  " << iter + 1 
      << "  SCF Iterations" << endl;
    this->fileio_->out << bannerEnd <<endl;
  }
};

