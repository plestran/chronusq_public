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
/***************************************
 * Compute Multipole Moments (COMPLEX) *
 ***************************************/
namespace ChronusQ {
template<>
void SingleSlater<dcomplex>::computeMultipole(){
  if(!this->haveDensity) this->formDensity();
  if(!this->aointegrals_->haveAOOneE) this->aointegrals_->computeAOOneE();
  if(!this->controls_->doDipole && !this->controls_->doQuadpole) return;

  int NB = this->nTCS_*this->nBasis_;
  int NBSq = NB*NB;
  int iBuf = 0;
  for(auto ixyz = 0; ixyz < 3; ixyz++){
    ConstRealMap mu(&this->aointegrals_->elecDipole_->storage()[iBuf],NB,NB);
    (*dipole_)(ixyz,0) = -this->densityA_->real().frobInner(mu);
    if(!this->isClosedShell && this->Ref_ != TCS) 
      (*dipole_)(ixyz,0) += -this->densityB_->real().frobInner(mu);
    iBuf += NBSq;
  }
  for(int iA = 0; iA < this->molecule_->nAtoms(); iA++)
    *this->dipole_ += elements[this->molecule_->index(iA)].atomicNumber *
          this->molecule_->cart()->col(iA);
  
  if(this->controls_->doQuadpole){
    iBuf = 0;
    for(auto jxyz = 0; jxyz < 3; jxyz++)
    for(auto ixyz = jxyz; ixyz < 3; ixyz++){
      ConstRealMap mu(&this->aointegrals_->elecQuadpole_->storage()[iBuf],NB,NB);
        (*quadpole_)(jxyz,ixyz) = -this->densityA_->real().frobInner(mu);
        if(!this->isClosedShell && this->Ref_ != TCS) 
          (*quadpole_)(jxyz,ixyz) += -this->densityB_->real().frobInner(mu);
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
      ConstRealMap mu(&this->aointegrals_->elecOctpole_->storage()[iBuf],NB,NB);
      (*octpole_)(kxyz,jxyz,ixyz) = -this->densityA_->real().frobInner(mu);
      if(!this->isClosedShell && this->Ref_ != TCS) 
        (*octpole_)(kxyz,jxyz,ixyz) += -this->densityB_->real().frobInner(mu);
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
//this->printMultipole();

}
template<>
void SingleSlater<dcomplex>::computeSExpect(){
  if(this->Ref_ == RHF){

    this->Sx_  = 0.0;
    this->Sy_  = 0.0;
    this->Sz_  = 0.0;
    this->Ssq_ = 0.0; 

  } else if(this->Ref_ == UHF) {

    this->Sx_  = 0.0;
    this->Sy_  = 0.0;
    this->Sz_  = 0.5*(this->nOccA_ - this->nOccB_);
    this->Ssq_ = this->Sz_ * (this->Sz_ + 1) + this->nOccB_;

    for(auto i = 0; i < this->nOccA_; i++)
    for(auto j = 0; j < this->nOccB_; j++){

      auto Sij = this->moA_->col(i).dot((*this->aointegrals_->overlap_) * this->moB_->col(j));
      this->Ssq_ -= std::real(Sij*std::conj(Sij));

    }

  } else if(this->Ref_ == CUHF) {

    this->Sx_  = 0.0;
    this->Sy_  = 0.0;
    this->Sz_  = 0.5*(this->nOccA_ - this->nOccB_);
    this->Ssq_ = this->Sz_ * (this->Sz_ + 1) + this->nOccB_;

  } else if(this->Ref_ == TCS) {
    
    this->Sx_  = 0.0;
    this->Sy_  = 0.0;
    this->Sz_  = 0.0;
    this->Ssq_  = 0.0;

    for(auto i = 0; i < this->nOccA_+this->nOccB_; i++) 
    for(auto j = 0; j < this->nOccA_+this->nOccB_; j++) {
      dcomplex SAA = 0.0;
      dcomplex SAB = 0.0;
      dcomplex SBB = 0.0;
      for(auto mu = 0; mu < this->nTCS_*this->nBasis_; mu += 2)
      for(auto nu = 0; nu < this->nTCS_*this->nBasis_; nu += 2){
        SAA += (*this->moA_)(mu,i) * 
               (*this->aointegrals_->overlap_)(mu,nu) * 
               (*this->moA_)(nu,j);

        SAB += (*this->moA_)(mu,i) * 
               (*this->aointegrals_->overlap_)(mu,nu) * 
               (*this->moA_)(nu+1,j);

        SBB += (*this->moA_)(mu+1,i) * 
               (*this->aointegrals_->overlap_)(mu,nu) * 
               (*this->moA_)(nu+1,j);
      }
      if( i == j ) {
        this->Sx_ += std::real(SAB);
        this->Sy_ -= std::imag(SAB);
        this->Sz_ += 0.5*std::real(SAA - SBB);
      }
      this->Ssq_ -= std::real(SAB*std::conj(SAB));
      this->Ssq_ -= 0.25*std::real((SAA-SBB)*std::conj(SAA-SBB));
    }
    this->Ssq_ += 0.75 * (this->nOccA_+this->nOccB_);
    this->Ssq_ += this->Sx_*this->Sx_;
    this->Ssq_ += this->Sy_*this->Sy_;
    this->Ssq_ += this->Sz_*this->Sz_;

  }
};
}; // namespace ChronusQ
