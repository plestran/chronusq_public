/* 
 *  This file is part of the Chronus Quantum (ChronusQ) software package
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

template<typename T>
void SingleSlater<T>::computeEnergy(){
/*
  if(getRank() == 0) {
    this->energyOneE = 
      this->template computeProperty<double,DENSITY_TYPE::TOTAL>(
          *this->aointegrals_->coreH_);
    if(this->nTCS_ == 1 && this->isClosedShell)
      this->energyTwoE = 
        0.5 * this->template computeProperty<double,DENSITY_TYPE::TOTAL>(
            this->PTA_->conjugate());
    else {
      this->energyTwoE = 
        0.5 * this->template computeProperty<double,DENSITY_TYPE::TOTAL>(
            this->PTScalar_->conjugate());
      this->energyTwoE +=
        0.5 * this->template computeProperty<double,DENSITY_TYPE::MZ>(
            this->PTMz_->conjugate());
      if(this->nTCS_ == 2) {
        this->energyTwoE +=
          0.5 * this->template computeProperty<double,DENSITY_TYPE::MX>(
              this->PTMx_->conjugate());
        this->energyTwoE +=
          0.5 * this->template computeProperty<double,DENSITY_TYPE::MY>(
              this->PTMy_->conjugate());
      }
      this->energyTwoE *= 0.5; // ??
    }
*/
  this->energyOneE = 
    this->template computeProperty<double,DENSITY_TYPE::TOTAL>(
        *this->aointegrals_->coreH_);
  
  if(this->aointegrals_->doX2C) {
    dcomplex SOPart(0);
    if(this->nTCS_ == 2 or !this->isClosedShell) {
      SOPart +=  this->template computeProperty<dcomplex,DENSITY_TYPE::MZ>(
            *this->aointegrals_->oneEmz_);
    }
    if(this->nTCS_ == 2){
      SOPart +=  this->template computeProperty<dcomplex,DENSITY_TYPE::MX>(
            *this->aointegrals_->oneEmx_);
      SOPart +=  this->template computeProperty<dcomplex,DENSITY_TYPE::MY>(
            *this->aointegrals_->oneEmy_);
    }
    SOPart *= math.ii;
    this->energyOneE -= std::real(SOPart);
  }

  this->energyTwoE = 
    0.5 * this->template computeProperty<double,DENSITY_TYPE::TOTAL>(
        this->PTScalar_->conjugate());
  if(this->nTCS_ == 2 or !this->isClosedShell) {
    this->energyTwoE += 
      0.5 * this->template computeProperty<double,DENSITY_TYPE::MZ>(
          this->PTMz_->conjugate());
  }
  if(this->nTCS_ == 2) {
    this->energyTwoE += 
      0.5 * this->template computeProperty<double,DENSITY_TYPE::MY>(
          this->PTMy_->conjugate());
    this->energyTwoE += 
      0.5 * this->template computeProperty<double,DENSITY_TYPE::MX>(
          this->PTMx_->conjugate());
  }
  this->energyTwoE *= 0.5;

      
    if(this->isDFT) this->energyTwoE += this->energyExc;

    // Add in the electric field component if they are non-zero
    std::array<double,3> null{{0,0,0}};
    if(this->elecField_ != null){
      this->computeMultipole();
      for(auto iXYZ = 0; iXYZ < 3; iXYZ++){
        this->energyOneE += this->elecField_[iXYZ] * this->elecDipole_[iXYZ];
      }
    }
 
    this->totalEnergy_= 
      this->energyOneE + this->energyTwoE + this->energyNuclei_;
#ifdef CQ_ENABLE_MPI
  MPI_Bcast(&this->totalEnergy_,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
  MPI_Bcast(&this->energyOneE,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
  MPI_Bcast(&this->energyTwoE,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
#endif
//this->gatherDensity();
};

template<typename T>
void SingleSlater<T>::computeMultipole(){
  if(!this->aointegrals_->haveAOOneE && getRank() == 0) 
    this->aointegrals_->computeAOOneE();
  if(this->maxMultipole_ < 1) return;

  if(getRank() == 0) {
    this->computeElecMultipoles(this->aointegrals_->elecDipoleSep_,
        this->aointegrals_->elecQuadpoleSep_,
        this->aointegrals_->elecOctpoleSep_);
    if(this->maxMultipole_ >= 1) {
      RealVecMap Dipole(&this->elecDipole_[0],3);
      for(int iA = 0; iA < this->molecule_->nAtoms(); iA++)
        Dipole += elements[this->molecule_->index(iA)].atomicNumber *
              this->molecule_->cart()->col(iA);
    }
    
    if(this->maxMultipole_ >= 2){
      RealMap QP(&this->elecQuadpole_[0][0],3,3);
      RealMap TQP(&this->elecTracelessQuadpole_[0][0],3,3);
      for(int iA = 0; iA < this->molecule_->nAtoms(); iA++){
        QP += 
          elements[this->molecule_->index(iA)].atomicNumber *
          this->molecule_->cart()->col(iA) * 
          this->molecule_->cart()->col(iA).transpose();
        TQP += 
          elements[this->molecule_->index(iA)].atomicNumber *
          this->molecule_->cart()->col(iA) * 
          this->molecule_->cart()->col(iA).transpose()
          - elements[this->molecule_->index(iA)].atomicNumber *
          RealMatrix::Identity(3,3) * 
          (this->molecule_->cart()->col(iA) * 
          this->molecule_->cart()->col(iA).transpose()).trace() /3;
      }
    }
 
    if(this->maxMultipole_ >= 3){
      for(auto iA = 0; iA < this->molecule_->nAtoms(); iA++)
      for(auto ixyz = 0; ixyz < 3; ixyz++)
      for(auto jxyz = 0; jxyz < 3; jxyz++)
      for(auto kxyz = 0; kxyz < 3; kxyz++)
        this->elecOctpole_[ixyz][jxyz][kxyz] += 
              elements[this->molecule_->index(iA)].atomicNumber *
              (*this->molecule_->cart())(ixyz,iA)*
              (*this->molecule_->cart())(jxyz,iA)*
              (*this->molecule_->cart())(kxyz,iA);
      
    }
  }

#ifdef CQ_ENABLE_MPI
  MPI_Bcast(&this->elecDipole_[0],3,MPI_DOUBLE,0,MPI_COMM_WORLD);
  MPI_Bcast(&this->elecQuadpole_[0][0],9,MPI_DOUBLE,0,MPI_COMM_WORLD);
  MPI_Bcast(&this->elecTracelessQuadpole_[0][0],9,MPI_DOUBLE,0,MPI_COMM_WORLD);
  MPI_Bcast(&this->elecOctpole_[0][0][0],27,MPI_DOUBLE,0,MPI_COMM_WORLD);
#endif

}

template<typename T>
void SingleSlater<T>::mullikenPop() {
  double charge;
  this->mullPop_.clear();
  (*this->NBSqScratch_) = 
    (*this->onePDMScalar_) * this->aointegrals_->overlap_->template cast<T>(); 

  for (auto iAtm = 0; iAtm < this->molecule_->nAtoms(); iAtm++) {
    auto iBfSt = this->basisset_->mapCen2Bf(iAtm)[0];
    auto iSize = this->basisset_->mapCen2Bf(iAtm)[1];
    charge  = elements[this->molecule_->index(iAtm)].atomicNumber;
    charge -= std::real(this->NBSqScratch_->
      block(iBfSt,iBfSt,iSize,iSize).trace());
    this->mullPop_.push_back(charge); 
  } 
}

template<typename T>
void SingleSlater<T>::loewdinPop() {
  double charge;
  this->lowPop_.clear();
  for (auto iAtm = 0; iAtm < this->molecule_->nAtoms(); iAtm++) {
    auto iBfSt = this->basisset_->mapCen2Bf(iAtm)[0];
    auto iSize = this->basisset_->mapCen2Bf(iAtm)[1];
    charge  = elements[this->molecule_->index(iAtm)].atomicNumber;
    charge -= std::real(this->onePDMOrthoScalar_->
      block(iBfSt,iBfSt,iSize,iSize).trace());
    this->lowPop_.push_back(charge); 
  } 
}
