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
  if(getRank() == 0) {
    if(!this->haveMO)
      CErr("No MO coefficients available to form one-particle density matrix!",
           this->fileio_->out);
 
    if(this->Ref_ == TCS){
      auto nOcc = this->nOccA_ + this->nOccB_;
      (*this->densityA_)= 
        this->moA_->block(0,0,this->nTCS_*this->nBasis_,nOcc)*
        this->moA_->block(0,0,this->nTCS_*this->nBasis_,nOcc).adjoint();
    } else {
      (*this->densityA_) = 
        this->moA_->block(0,0,this->nBasis_,this->nOccA_)*
        this->moA_->block(0,0,this->nBasis_,this->nOccA_).adjoint();
      // D(a) is actually total D for RHF
      if(this->isClosedShell) (*this->densityA_) *= math.two;
      else {
        (*this->densityB_) = 
          this->moB_->block(0,0,this->nBasis_,this->nOccB_)*
          this->moB_->block(0,0,this->nBasis_,this->nOccB_).adjoint();
      }
    }
    if(this->printLevel_ >= 2) this->printDensity();
  }
  this->haveDensity = true;
  this->mpiBCastDensity();

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
void SingleSlater<T>::mpiBCastDensity(){
#ifdef CQ_ENABLE_MPI
  auto dataType = MPI_DOUBLE;
  if(typeid(T).hash_code() == typeid(dcomplex).hash_code())
    dataType = MPI_C_DOUBLE_COMPLEX;

  MPI_Bcast(this->densityA_->data(),
    this->nTCS_*this->nTCS_*this->nBasis_*this->nBasis_,dataType,0,
    MPI_COMM_WORLD);
  if(!this->isClosedShell && this->Ref_ != TCS)
    MPI_Bcast(this->densityB_->data(),
      this->nTCS_*this->nTCS_*this->nBasis_*this->nBasis_,dataType,0,
      MPI_COMM_WORLD);
#endif
  ;
};
