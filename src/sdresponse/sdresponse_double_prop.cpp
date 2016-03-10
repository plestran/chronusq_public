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
void SDResponse<double>::formTransDipole(){
   auto NBSq = this->nTCS_*this->nTCS_*this->nBasis_*this->nBasis_;
   if(this->iMeth_ == CIS || this->iMeth_ == RPA || this->iMeth_ == STAB){
     if(this->Ref_ == SingleSlater<double>::TCS) {
       RealMatrix TAO(this->nTCS_*this->nBasis_,this->nTCS_*this->nBasis_);
       for(auto iSt = 0; iSt < this->nSek_; iSt++){
         RealVecMap TMOV(this->transDen_->data()+iSt*this->nSingleDim_,this->nSingleDim_);
         this->formAOTDen(TMOV,TAO,TAO);
         for(auto iXYZ = 0, iOff = 0; iXYZ < 3; iXYZ++, iOff += NBSq){
           RealMap dipole(&this->elecDipole_->storage()[iOff],this->nTCS_*this->nBasis_,this->nTCS_*this->nBasis_);
           (*this->transDipole_)(0,iSt+1,iXYZ) = TAO.frobInner(dipole);
         }
       }
     } else {
       RealMatrix TAOA(this->nBasis_,this->nBasis_);
       RealMatrix TAOB(this->nBasis_,this->nBasis_);
       for(auto iSt = 0; iSt < this->nSek_; iSt++){
         RealVecMap TMOV(this->transDen_->data()+iSt*this->nSingleDim_,this->nSingleDim_);
         this->formAOTDen(TMOV,TAOA,TAOB);
         for(auto iXYZ = 0, iOff = 0; iXYZ < 3; iXYZ++, iOff += NBSq){
           RealMap dipole(&this->elecDipole_->storage()[iOff],this->nBasis_,this->nBasis_);
           (*this->transDipole_)(0,iSt+1,iXYZ) = (TAOA+TAOB).frobInner(dipole);
         }
       }
     }
   } else if(this->iMeth_ == PPATDA) {
     if(this->Ref_ == SingleSlater<double>::TCS) {
       RealMatrix DMO(this->nTCS_*this->nBasis_,this->nTCS_*this->nBasis_);
       for(auto iXYZ = 0, iOff = 0; iXYZ < 3; iXYZ++, iOff += NBSq){
         RealMap dipole(&this->elecDipole_->storage()[iOff],this->nTCS_*this->nBasis_,this->nTCS_*this->nBasis_);
         cout << "HERE" << endl;
         DMO = this->singleSlater_->moA()->adjoint() * dipole * (*this->singleSlater_->moA());
         cout << "HERE" << endl;
         double occTrace = 0;
         for(auto i = 0; i < this->nO_; i++) occTrace += DMO(i,i);
         for(auto iSt = 0; iSt < this->nSek_; iSt++)
         for(auto jSt = 0; jSt < iSt;         jSt++){
           for(auto a = 0, ab = 0; a < this->nV_; a++)
           for(auto b = 0; b <  a ; b++, ab++    )
           for(auto c = 0, cd = 0; c < this->nV_; c++)
           for(auto d = 0; d <  c ; d++, cd++    ){
             auto prod = (*this->transDen_)(ab,iSt) * (*this->transDen_)(cd,jSt);
             if(a == c) (*this->transDipole_)(iSt,jSt,iXYZ) += prod*DMO(this->nO_+d,this->nO_+b);
             if(a == d) (*this->transDipole_)(iSt,jSt,iXYZ) -= prod*DMO(this->nO_+c,this->nO_+b);
             if(b == c) (*this->transDipole_)(iSt,jSt,iXYZ) -= prod*DMO(this->nO_+d,this->nO_+a);
             if(b == d) (*this->transDipole_)(iSt,jSt,iXYZ) += prod*DMO(this->nO_+c,this->nO_+a);
             if(a == c || b == d) (*this->transDipole_)(iSt,jSt,iXYZ) += occTrace*prod;
             if(a == d || b == c) (*this->transDipole_)(iSt,jSt,iXYZ) -= occTrace*prod;
           }
           (*this->transDipole_)(jSt,iSt,iXYZ) = (*this->transDipole_)(iSt,jSt,iXYZ);
         }
       }
     } else CErr("Transition Dipole for non-TCS ppTDA NYI",this->fileio_->out); 
   } else CErr("Transition Dipole not defined to PSCF Method",this->fileio_->out);
} //formTransDipole

template<>
void SDResponse<double>::formOscStrength(){
  this->oscStrength_->setZero();
  if(this->iMeth_ == PPATDA){
    for(auto iSt  = 0; iSt  < this->nSek_; iSt++ )
    for(auto jSt  = 0; jSt  < iSt        ; jSt++ ){
      for(auto iXYZ = 0; iXYZ < 3;           iXYZ++){
        (*this->oscStrength_)(iSt,jSt) +=
          (2.0/3.0)*((*this->omega_)(iSt) - (*this->omega_)(jSt))*
          std::pow((*this->transDipole_)(iSt,jSt,iXYZ),2.0);
      }
      (*this->oscStrength_)(jSt,iSt) = (*this->oscStrength_)(iSt,jSt);
    }
  } else {
    for(auto iSt  = 0; iSt  < this->nSek_; iSt++ )
    for(auto iXYZ = 0; iXYZ < 3;           iXYZ++){
      (*this->oscStrength_)(0,iSt+1) +=
        (2.0/3.0)*(*this->omega_)(iSt)*
        std::pow((*this->transDipole_)(0,iSt+1,iXYZ),2.0);
    }
  }
} //formOscStrength
} // namespace ChronusQ
