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

template<typename T>
void SingleSlater<T>::computeEnergy(){
  if(getRank() == 0) {
    this->energyOneE = 
      this->template computeProperty<double,DENSITY_TYPE::TOTAL>(
          *this->aointegrals_->oneE_);
    if(!this->isClosedShell && this->Ref_ != TCS)
      this->energyTwoE = 0.5 * (
        this->template computeProperty<double,DENSITY_TYPE::ALPHA>(this->PTA_->conjugate()) 
        +this->template computeProperty<double,DENSITY_TYPE::BETA>(this->PTB_->conjugate())
      );
    else
      this->energyTwoE = 0.5 * (
        this->template computeProperty<double,DENSITY_TYPE::TOTAL>(this->PTA_->conjugate()) 
      ); 
      
    if(this->isDFT) this->energyTwoE += this->totalEx + this->totalEcorr;

    // Add in the electric field component if they are non-zero
    std::array<double,3> null{{0,0,0}};
    if(this->elecField_ != null){
    //int NB = this->nTCS_*this->nBasis_;
    //int NBSq = NB*NB;
    //int iBuf = 0;
      auto exptdipole = this-> template computeProperty<double,
           DENSITY_TYPE::TOTAL>(this->aointegrals_->elecDipoleSep_);
      for(auto iXYZ = 0; iXYZ < 3; iXYZ++){
      //RealMap mu(&this->aointegrals_->elecDipole_->storage()[iBuf],NB,NB);
      //this->energyOneE += this->elecField_[iXYZ] *
      //  this->template computeProperty<double,DENSITY_TYPE::TOTAL>(mu);
      //iBuf += NBSq;
        this->energyOneE += this->elecField_[iXYZ] * exptdipole[iXYZ];
      }
    }
 
 
    this->totalEnergy= this->energyOneE + this->energyTwoE + this->energyNuclei;
  }
#ifdef CQ_ENABLE_MPI
  MPI_Bcast(&this->totalEnergy,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
  MPI_Bcast(&this->energyOneE,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
  MPI_Bcast(&this->energyTwoE,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
#endif
//this->printEnergy();
};

template<typename T>
void SingleSlater<T>::computeMultipole(){
  if(!this->haveDensity) this->formDensity();
  if(!this->aointegrals_->haveAOOneE && getRank() == 0) 
    this->aointegrals_->computeAOOneE();
  if(this->maxMultipole_ < 1) return;

  if(getRank() == 0) {
    int NB = this->nTCS_*this->nBasis_;
    int NBSq = NB*NB;
    int iBuf = 0;
 
    if(this->maxMultipole_ >= 1) {
      auto exptdipole = this-> template computeProperty<double,
           DENSITY_TYPE::TOTAL>(this->aointegrals_->elecDipoleSep_);

      for(int iXYZ = 0; iXYZ < 3; iXYZ++)
        (*this->dipole_)(iXYZ,0) = - exptdipole[iXYZ];

      for(int iA = 0; iA < this->molecule_->nAtoms(); iA++)
        *this->dipole_ += elements[this->molecule_->index(iA)].atomicNumber *
              this->molecule_->cart()->col(iA);
    }
    
    if(this->maxMultipole_ >= 2){
      auto exptquadpole = this-> template computeProperty<double,
           DENSITY_TYPE::TOTAL>(this->aointegrals_->elecQuadpoleSep_);

      for(auto jxyz = 0, iX = 0; jxyz < 3; jxyz++)
      for(auto ixyz = jxyz; ixyz < 3; ixyz++, iX++)
        (*this->quadpole_)(jxyz,ixyz) = - exptquadpole[iX];

      (*this->quadpole_) = this->quadpole_->template selfadjointView<Upper>();

      for(int iA = 0; iA < this->molecule_->nAtoms(); iA++)
        *this->quadpole_ += elements[this->molecule_->index(iA)].atomicNumber *
              this->molecule_->cart()->col(iA) * 
              this->molecule_->cart()->col(iA).transpose();
 
      (*this->tracelessQuadpole_) = (*this->quadpole_) - 
        RealMatrix::Identity(3,3)*this->quadpole_->trace()/3;
    }
 
    if(this->maxMultipole_ >= 3){
      auto exptOctpole = this-> template computeProperty<double,
           DENSITY_TYPE::TOTAL>(this->aointegrals_->elecOctpoleSep_);

      for(auto kxyz = 0,iX = 0;    kxyz < 3; kxyz++)
      for(auto jxyz = kxyz; jxyz < 3; jxyz++)
      for(auto ixyz = jxyz; ixyz < 3; ixyz++, iX++)
        (*this->octpole_)(kxyz,jxyz,ixyz) = - exptOctpole[iX]; 

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
        (*this->octpole_)(ixyz,jxyz,kxyz) += 
              elements[this->molecule_->index(iA)].atomicNumber *
              (*this->molecule_->cart())(ixyz,iA)*
              (*this->molecule_->cart())(jxyz,iA)*
              (*this->molecule_->cart())(kxyz,iA);
      
    }
  }
#ifdef CQ_ENABLE_MPI
  MPI_Bcast(this->dipole_->data(),3,MPI_DOUBLE,0,MPI_COMM_WORLD);
  MPI_Bcast(this->quadpole_->data(),9,MPI_DOUBLE,0,MPI_COMM_WORLD);
  MPI_Bcast(this->tracelessQuadpole_->data(),9,MPI_DOUBLE,0,MPI_COMM_WORLD);
  MPI_Bcast(this->octpole_->data(),27,MPI_DOUBLE,0,MPI_COMM_WORLD);
#endif

}

template<typename T>
void SingleSlater<T>::mullikenPop() {
  double charge;
  this->mullPop_.clear();
  RealMatrix PS = (*this->densityA_).real() * (*this->aointegrals_->overlap_); 
  if(!this->isClosedShell && this->Ref_ != TCS){ 
    PS += (*this->densityB_).real() * (*this->aointegrals_->overlap_);
  }
  for (auto iAtm = 0; iAtm < this->molecule_->nAtoms(); iAtm++) {
    auto iBfSt = this->basisset_->mapCen2Bf(iAtm)[0];
    auto iSize = this->basisset_->mapCen2Bf(iAtm)[1];
    charge  = elements[this->molecule_->index(iAtm)].atomicNumber;
    charge -= PS.block(iBfSt,iBfSt,iSize,iSize).trace();
    this->mullPop_.push_back(charge); 
  } 
}

template<typename T>
void SingleSlater<T>::computeSExpect(){
  if(getRank() == 0){
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
 
        dcomplex Sij = 
          this->moA_->col(i).dot(
            (*this->aointegrals_->overlap_) * this->moB_->col(j)
          );

        this->Ssq_ -= std::real(Sij*std::conj(Sij));
 
      }
 
    } else if(this->Ref_ == CUHF) {
 
      this->Sx_  = 0.0;
      this->Sy_  = 0.0;
      this->Sz_  = 0.5*(this->nOccA_ - this->nOccB_);
      this->Ssq_ = this->Sz_ * (this->Sz_ + 1);
 
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
  }
#ifdef CQ_ENABLE_MPI
  MPI_Bcast(&this->Sx_,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
  MPI_Bcast(&this->Sy_,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
  MPI_Bcast(&this->Sz_,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
  MPI_Bcast(&this->Ssq_,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
#endif
};

