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
#include <singleslater.h>
/************************************
 * Compute Multipole Moments (REAL) *
 ************************************/
namespace ChronusQ {
template<>
void SingleSlater<double>::computeMultipole(){
  if(!this->haveDensity) this->formDensity();
  if(!this->aointegrals_->haveAOOneE) this->aointegrals_->computeAOOneE();
  if(!this->controls_->doDipole && !this->controls_->doQuadpole) return;

  int NB = this->nBasis_;
  int NBSq = NB*NB;
  int iBuf = 0;
  for(auto ixyz = 0; ixyz < 3; ixyz++){
    ConstRealMap mu(&this->aointegrals_->elecDipole_->storage()[iBuf],NB,NB);
    (*dipole_)(ixyz,0) = -this->densityA_->frobInner(mu);
    iBuf += NBSq;
  }
  for(int iA = 0; iA < this->molecule_->nAtoms(); iA++)
    *this->dipole_ += elements[this->molecule_->index(iA)].atomicNumber *
          this->molecule_->cart()->col(iA);
  
  if(this->controls_->doQuadpole){
    iBuf = 0;
    for(auto jxyz = 0; jxyz < 3; jxyz++)
    for(auto ixyz = jxyz; ixyz < 3; ixyz++){
      ConstRealMap 
        mu(&this->aointegrals_->elecQuadpole_->storage()[iBuf],NB,NB);
      (*quadpole_)(jxyz,ixyz) = -this->densityA_->frobInner(mu);
      iBuf += NBSq;
    }
    *this->quadpole_ = this->quadpole_->selfadjointView<Upper>();
    for(int iA = 0; iA < this->molecule_->nAtoms(); iA++)
      *this->quadpole_ += elements[this->molecule_->index(iA)].atomicNumber *
            this->molecule_->cart()->col(iA) * 
            this->molecule_->cart()->col(iA).transpose();

    (*this->tracelessQuadpole_) = (*this->quadpole_) - RealMatrix::Identity(3,3)*this->quadpole_->trace()/3;
  }

  if(this->controls_->doOctpole){
    iBuf = 0;
    for(auto kxyz = 0;    kxyz < 3; kxyz++)
    for(auto jxyz = kxyz; jxyz < 3; jxyz++)
    for(auto ixyz = jxyz; ixyz < 3; ixyz++){
      ConstRealMap 
        mu(&this->aointegrals_->elecOctpole_->storage()[iBuf],NB,NB);
      (*octpole_)(kxyz,jxyz,ixyz) = -this->densityA_->frobInner(mu);
      iBuf += NBSq;
    }
    for(auto kxyz = 0;    kxyz < 3; kxyz++)
    for(auto jxyz = kxyz; jxyz < 3; jxyz++)
    for(auto ixyz = jxyz; ixyz < 3; ixyz++){
      (*octpole_)(ixyz,jxyz,kxyz) = (*octpole_)(kxyz,jxyz,ixyz);
      (*octpole_)(ixyz,kxyz,jxyz) = (*octpole_)(kxyz,jxyz,ixyz);
      (*octpole_)(jxyz,ixyz,kxyz) = (*octpole_)(kxyz,jxyz,ixyz);
      (*octpole_)(jxyz,kxyz,ixyz) = (*octpole_)(kxyz,jxyz,ixyz);
      (*octpole_)(kxyz,ixyz,jxyz) = (*octpole_)(kxyz,jxyz,ixyz);
    }

    for(auto iA = 0; iA < this->molecule_->nAtoms(); iA++)
    for(auto ixyz = 0; ixyz < 3; ixyz++)
    for(auto jxyz = 0; jxyz < 3; jxyz++)
    for(auto kxyz = 0; kxyz < 3; kxyz++)
      (*this->octpole_)(ixyz,jxyz,kxyz) += elements[this->molecule_->index(iA)].atomicNumber *
            (*this->molecule_->cart())(ixyz,iA)*
            (*this->molecule_->cart())(jxyz,iA)*
            (*this->molecule_->cart())(kxyz,iA);
    
  }
  this->printMultipole();

}
}; // namespace ChronusQ
