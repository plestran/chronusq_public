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
void RealTime<T>::iniRealTime(Molecule * molecule, BasisSet *basisset, FileIO *fileio, Controls *controls, AOIntegrals *aointegrals,SingleSlater<T> *groundState) {

  this->fileio_         = fileio;
  this->controls_       = controls;
  this->aointegrals_	= aointegrals;
  this->groundState_   	= groundState;

  this->ssPropagator_	= std::unique_ptr<SingleSlater<dcomplex>>(new SingleSlater<dcomplex>(groundState));

  this->nBasis_ = basisset->nBasis();
  this->RHF_	= this->groundState_->isClosedShell;
  this->nOccA_ 	= this->groundState_->nOccA();
  this->nOccB_ 	= this->groundState_->nOccB();

  this->frozenNuc_	= true;
  this->maxSteps_	= this->controls_->rtMaxSteps;
  this->stepSize_	= this->controls_->rtTimeStep;
  this->typeOrtho_	= this->controls_->rtTypeOrtho; 
  this->initDensity_	= this->controls_->rtInitDensity;
  this->swapMOA_	= this->controls_->rtSwapMOA;
  this->swapMOB_	= this->controls_->rtSwapMOB;
  this->methFormU_	= this->controls_->rtMethFormU;

  this->fileio_->out<<"\nReal-time TDHF: "<<endl;
  this->fileio_->out<<std::right<<std::setw(20)<<"Number of steps = "<<std::setw(15)<<this->maxSteps_<<std::setw(5)<<endl;
  this->fileio_->out<<std::right<<std::setw(20)<<"Step size = "<<std::setw(15)<<this->stepSize_<<std::setw(5)<<" a.u. "<<endl;

  this->oTrans1_ = std::unique_ptr<ComplexMatrix>(new ComplexMatrix(this->nBasis_,this->nBasis_));
  this->oTrans2_ = std::unique_ptr<ComplexMatrix>(new ComplexMatrix(this->nBasis_,this->nBasis_));
  this->POA_  	 = std::unique_ptr<ComplexMatrix>(new ComplexMatrix(this->nBasis_,this->nBasis_));
  this->POAsav_  = std::unique_ptr<ComplexMatrix>(new ComplexMatrix(this->nBasis_,this->nBasis_));
  this->POB_  	 = std::unique_ptr<ComplexMatrix>(new ComplexMatrix(this->nBasis_,this->nBasis_));
  this->POBsav_  = std::unique_ptr<ComplexMatrix>(new ComplexMatrix(this->nBasis_,this->nBasis_));
  this->FOA_ 	 = std::unique_ptr<ComplexMatrix>(new ComplexMatrix(this->nBasis_,this->nBasis_));
  this->FOB_ 	 = std::unique_ptr<ComplexMatrix>(new ComplexMatrix(this->nBasis_,this->nBasis_));
  this->initMOA_ = std::unique_ptr<ComplexMatrix>(new ComplexMatrix(this->nBasis_,this->nBasis_));
  this->initMOB_ = std::unique_ptr<ComplexMatrix>(new ComplexMatrix(this->nBasis_,this->nBasis_));
  this->uTransA_ = std::unique_ptr<ComplexMatrix>(new ComplexMatrix(this->nBasis_,this->nBasis_));
  this->uTransB_ = std::unique_ptr<ComplexMatrix>(new ComplexMatrix(this->nBasis_,this->nBasis_));
  this->scratch_ = std::unique_ptr<ComplexMatrix>(new ComplexMatrix(this->nBasis_,this->nBasis_));
};

