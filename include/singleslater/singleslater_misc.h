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
/***********************
 * Form Density Matrix *
 ***********************/
template<typename T>
void SingleSlater<T>::formDensity(){
  if(getRank() == 0) {
    if(!this->haveMO)
      CErr("No MO coefficients available to form one-particle density matrix!",
           this->fileio_->out);
 
    cout << "HERE 9" << endl;
    if(this->Ref_ == TCS){
      auto nOcc = this->nOccA_ + this->nOccB_;
      (*this->onePDMA_) = 
        this->moA_->block(0,0,this->nTCS_*this->nBasis_,nOcc)*
        this->moA_->block(0,0,this->nTCS_*this->nBasis_,nOcc).adjoint();
    } else {
    cout << "HERE 9" << endl;
      (*this->onePDMA_) = 
        this->moA_->block(0,0,this->nBasis_,this->nOccA_)*
        this->moA_->block(0,0,this->nBasis_,this->nOccA_).adjoint();
      // D(a) is actually total D for RHF
      if(this->isClosedShell) (*this->onePDMA_) *= math.two;
      else {
    cout << "HERE 9" << endl;
        (*this->onePDMB_) = 
          this->moB_->block(0,0,this->nBasis_,this->nOccB_)*
          this->moB_->block(0,0,this->nBasis_,this->nOccB_).adjoint();
      }
    }
    if(this->printLevel_ >= 2) this->printDensity();
  }
  this->haveDensity = true;
  this->mpiBCastDensity();
}
