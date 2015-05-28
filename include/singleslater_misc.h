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
/***********************
 * Form Density Matrix *
 ***********************/
template<typename T>
void SingleSlater<T>::formDensity(){
  if(!this->haveMO)
    CErr("No MO coefficients available to form one-particle density matrix!",
         this->fileio_->out);
  *this->densityA_ = this->moA_->block(0,0,this->nBasis_,this->nOccA_)*
                   this->moA_->block(0,0,this->nBasis_,this->nOccA_).adjoint();
  if(this->RHF_) *this->densityA_ *= math.two;
  else {
    *this->densityB_ = this->moB_->block(0,0,this->nBasis_,this->nOccB_)*
                   this->moB_->block(0,0,this->nBasis_,this->nOccB_).adjoint();
  }
  if(this->controls_->printLevel>=2) {
    prettyPrint(this->fileio_->out,(*this->densityA_),"Alpha Density");
    if(!this->RHF_) prettyPrint(this->fileio_->out,(*this->densityB_),"Beta Density");
  };
  this->haveDensity = true;
}
/************************
 * Compute Total Energy *
 ************************/
template<typename T>
void SingleSlater<T>::computeEnergy(){
  this->energyOneE = (*this->aointegrals_->oneE_).frobInner(this->densityA_->conjugate());
#ifndef USE_LIBINT
  this->energyTwoE = ((*this->coulombA_)-(*this->exchangeA_)).frobInner(this->densityA->conjugate());
#else
  this->energyTwoE = 0.5*(*this->PTA_).frobInner(this->densityA_->conjugate());
#endif
  this->totalEnergy= this->energyOneE + this->energyTwoE + this->energyNuclei;
  this->printEnergy();
};

