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
#include <sdresponse.h>
using ChronusQ::SDResponse;

namespace ChronusQ {
template<>
void SDResponse<dcomplex>::formTransDipole(){
   ComplexMatrix TAOA(this->nBasis_,this->nBasis_);
   ComplexMatrix TAOB(this->nBasis_,this->nBasis_);
   auto NBSq = this->nBasis_*this->nBasis_;
   for(auto iSt = 0; iSt < this->nSek_; iSt++){
     ComplexVecMap TMOV(this->transDen_->data()+iSt*this->nSingleDim_,this->nSingleDim_);
     this->formAOTDen(TMOV,TAOA,TAOB);
     for(auto iXYZ = 0, iOff = 0; iXYZ < 3; iXYZ++, iOff += NBSq){
       RealMap dipole(&this->elecDipole_->storage()[iOff],this->nBasis_,this->nBasis_);
       (*this->transDipole_)(0,iSt+1,iXYZ) = (TAOA+TAOB).real().frobInner(dipole);
     }
   }
} //formTransDipole

template<>
void SDResponse<dcomplex>::formOscStrength(){
  this->oscStrength_->setZero();
  for(auto iSt  = 0; iSt  < this->nSek_; iSt++ )
  for(auto iXYZ = 0; iXYZ < 3;           iXYZ++){
    (*this->oscStrength_)(0,iSt+1) +=
      (2.0/3.0)*(*this->omega_)(iSt)*
      std::pow((*this->transDipole_)(0,iSt+1,iXYZ),2.0);
  }
} //formOscStrength
} // namespace ChronusQ
