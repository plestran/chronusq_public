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

  if(this->Ref_ == TCS){
    auto nOcc = this->nOccA_ + this->nOccB_;
    (*this->densityA_)= 
      this->moA_->block(0,0,this->nTCS_*this->nBasis_,nOcc)*
      this->moA_->block(0,0,this->nTCS_*this->nBasis_,nOcc).adjoint();
/*
    for(auto i = 0; i < this->nTCS_*this->nBasis_; i+=2)
    for(auto j = 0; j < this->nTCS_*this->nBasis_; j+=2){
      (*this->densityA_)(i,j)     = 0.0;
      (*this->densityA_)(i+1,j+1) = 0.0;
      (*this->densityA_)(i+1,j)   = 0.0;
      (*this->densityA_)(i,j+1)   = 0.0;
      for(auto k = 0; k < nOcc; k++){
        (*this->densityA_)(i,j)     += (*this->moA_)(i,k)   * (*this->moA_)(j,k);
        (*this->densityA_)(i+1,j+1) += (*this->moA_)(i+1,k) * (*this->moA_)(j+1,k);
        (*this->densityA_)(i+1,j)   += (*this->moA_)(i+1,k) * (*this->moA_)(j,k);
        (*this->densityA_)(i,j+1)   += (*this->moA_)(i,k)   * (*this->moA_)(j+1,k);
      }
    }
*/
  } else {
    (*this->densityA_) = 
      this->moA_->block(0,0,this->nBasis_,this->nOccA_)*
      this->moA_->block(0,0,this->nBasis_,this->nOccA_).adjoint();
    if(this->isClosedShell) (*this->densityA_) *= math.two;// D(a) is actually total D for RHF
    else {
      (*this->densityB_) = 
        this->moB_->block(0,0,this->nBasis_,this->nOccB_)*
        this->moB_->block(0,0,this->nBasis_,this->nOccB_).adjoint();
    }
  }
  if(this->controls_->printLevel>=2) this->printDensity();
  this->haveDensity = true;
}
