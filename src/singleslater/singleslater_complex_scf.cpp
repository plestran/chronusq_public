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
#include <singleslater.h>
using ChronusQ::FileIO;
using ChronusQ::SingleSlater;

//----------------------------------------//
// do the SCF                             //
// Sajan                                  //
//----------------------------------------//
namespace ChronusQ {
/*
template<>
void SingleSlater<dcomplex>::formNO(){
  int INFO;
  char JOBZ = 'V';
  char UPLO = 'L';

  ComplexMap P(this->PNOMem_,this->nBasis_,this->nBasis_);

  P.real() = 0.5 * (*this->aointegrals_->ortho2_) * this->onePDMA_->real() * 
    (*this->aointegrals_->ortho2_);
  P.imag() = 0.5 * (*this->aointegrals_->ortho2_) * this->onePDMA_->imag() * 
    (*this->aointegrals_->ortho2_);

  if(!this->isClosedShell){
    P.real() += 0.5 * (*this->aointegrals_->ortho2_) * this->onePDMB_->real() *
      (*this->aointegrals_->ortho2_);
    P.imag() += 0.5 * (*this->aointegrals_->ortho2_) * this->onePDMB_->imag() *
      (*this->aointegrals_->ortho2_);
  }

  int LWORK  = 5*this->nTCS_*this->nBasis_;
  int LRWORK = 3*this->nTCS_*this->nBasis_;
  dcomplex *WORK  = this->memManager_->malloc<dcomplex>(LWORK);
  double   *RWORK = this->memManager_->malloc<double>(LRWORK);

  zheev_(&JOBZ,&UPLO,&this->nBasis_,this->PNOMem_,&this->nBasis_,
      this->occNumMem_,WORK,&LWORK,RWORK,&INFO);

  this->memManager_->free(WORK,LWORK);
  this->memManager_->free(RWORK,LRWORK);

  if(INFO != 0) CErr("ZHEEV Failed in FormNO",this->fileio_->out);
//P.transposeInPlace(); //bc row major

  // Swap Ordering
  for(auto i = 0; i < this->nBasis_/2; i++) 
    P.col(i).swap(P.col(this->nBasis_ - i- 1));

}
*/
template<>
void SingleSlater<dcomplex>::formNO(){
  int INFO;
  char JOBZ = 'V';
  char UPLO = 'L';

  ComplexMap P(this->PNOMem_,this->nBasis_,this->nBasis_);

  this->aointegrals_->Ortho2Trans(*this->onePDMScalar_,P);
  P *= 0.5;

  int LWORK  = 5*this->nTCS_*this->nBasis_;
  int LRWORK = 3*this->nTCS_*this->nBasis_;
  dcomplex *WORK  = this->memManager_->malloc<dcomplex>(LWORK);
  double   *RWORK = this->memManager_->malloc<double>(LRWORK);

  zheev_(&JOBZ,&UPLO,&this->nBasis_,this->PNOMem_,&this->nBasis_,
      this->occNumMem_,WORK,&LWORK,RWORK,&INFO);

  this->memManager_->free(WORK,LWORK);
  this->memManager_->free(RWORK,LRWORK);

  if(INFO != 0) CErr("ZHEEV Failed in FormNO",this->fileio_->out);

  // Swap Ordering
  for(auto i = 0; i < this->nBasis_/2; i++) 
    P.col(i).swap(P.col(this->nBasis_ - i - 1));

}


template<>
void SingleSlater<dcomplex>::mixOrbitalsComplex(){
  this->fileio_->out << "** Mixing HOMO and LUMO for Complex Guess **" << endl;
//auto nO = this->nAE_ + this->nBE_;
  if (this->Ref_==TCS) {
    auto HOMO = this->moA_->col(this->nO_-1);
    auto LUMO = this->moA_->col(this->nO_);
    this->moA_->col(this->nO_-1) = std::sqrt(0.5) * (HOMO + math.ii*LUMO);
    this->moA_->col(this->nO_)   = std::sqrt(0.5) * (HOMO - math.ii*LUMO);
  } else {
    auto HOMO = this->moA_->col(this->nOA_-1);
    auto LUMO = this->moA_->col(this->nOA_);
    this->moA_->col(this->nOA_-1) = std::sqrt(0.5) * (HOMO + math.ii*LUMO);
    this->moA_->col(this->nOA_)   = std::sqrt(0.5) * (HOMO - math.ii*LUMO);
  }
}

template<>
void SingleSlater<dcomplex>::diagFock2(){
  int INFO;
  char JOBZ = 'V';
  char UPLO = 'U';
  auto NTCSxNBASIS = this->nTCS_*this->nBasis_;
  
  int LWORK  = 5*NTCSxNBASIS;
  int LRWORK = 3*NTCSxNBASIS;
  dcomplex *WORK  = this->memManager_->malloc<dcomplex>(LWORK);
  double   *RWORK = this->memManager_->malloc<double>(LRWORK);
/*
  zheev_(&JOBZ,&UPLO,&NTCSxNBASIS,this->fockOrthoA_->data(),&NTCSxNBASIS,
      this->epsA_->data(),WORK,&LWORK,RWORK,&INFO);
  if(INFO != 0) CErr("ZHEEV Failed Fock Alpha",this->fileio_->out);
  (*this->moA_) = (*this->fockOrthoA_);

  if(this->nTCS_ == 1 && !this->isClosedShell){
    zheev_(&JOBZ,&UPLO,&NTCSxNBASIS,this->fockOrthoB_->data(),&NTCSxNBASIS,
        this->epsB_->data(),WORK,&LWORK,RWORK,&INFO);
    if(INFO != 0) CErr("ZHEEV Failed Fock Beta",this->fileio_->out);
    (*this->moB_) = (*this->fockOrthoB_);
  }
*/
  // ** WARNING **
  // this assumes that the orthonormal fock has been copied into
  // the MO coefficients storage
  zheev_(&JOBZ,&UPLO,&NTCSxNBASIS,this->moA_->data(),&NTCSxNBASIS,
      this->epsA_->data(),WORK,&LWORK,RWORK,&INFO);
  if(INFO != 0) CErr("ZHEEV Failed Fock Alpha",this->fileio_->out);

  if(this->nTCS_ == 1 && !this->isClosedShell){
    zheev_(&JOBZ,&UPLO,&NTCSxNBASIS,this->moB_->data(),&NTCSxNBASIS,
        this->epsB_->data(),WORK,&LWORK,RWORK,&INFO);
    if(INFO != 0) CErr("ZHEEV Failed Fock Beta",this->fileio_->out);
  }
  this->memManager_->free(WORK,LWORK);
  this->memManager_->free(RWORK,LRWORK);
};


/*
template<>
void SingleSlater<dcomplex>::fockCUHF() {
  ComplexMap P(this->PNOMem_,this->nBasis_,this->nBasis_);
  ComplexMap DelF(this->delFMem_,this->nBasis_,this->nBasis_);
  ComplexMap Lambda(this->lambdaMem_,this->nBasis_,this->nBasis_);

  int activeSpace  = this->molecule_->multip() - 1;
  int coreSpace    = (this->molecule_->nTotalE() - activeSpace) / 2;
  int virtualSpace = this->nBasis_ - coreSpace - activeSpace;

  // DelF = X * (F(A) - F(B)) * X
  this->NBSqScratch_->real() = 0.5 * (*this->aointegrals_->ortho1_) *
    this->fockMz_->real();
  this->NBSqScratch_->imag() = 0.5 * (*this->aointegrals_->ortho1_) *
    this->fockMz_->imag();
  DelF.real() = this->NBSqScratch_->real() * (*this->aointegrals_->ortho1_);
  DelF.imag() = this->NBSqScratch_->imag() * (*this->aointegrals_->ortho1_);

  // DelF = C(NO)^\dagger * DelF * C(NO) (Natural Orbitals)
  (*this->NBSqScratch_) = P.adjoint() * DelF;
  DelF = (*this->NBSqScratch_) * P;

  Lambda.setZero();
  for(auto i = activeSpace + coreSpace; i < this->nBasis_; i++)
  for(auto j = 0                      ; j < coreSpace    ; j++){
    Lambda(i,j) = -DelF(i,j);
    Lambda(j,i) = -DelF(j,i);
  }

  (*this->NBSqScratch_) = P * Lambda;
  Lambda = (*this->NBSqScratch_) * P.transpose();

  this->NBSqScratch_->real() = (*this->aointegrals_->ortho2_) * Lambda.real();
  this->NBSqScratch_->imag() = (*this->aointegrals_->ortho2_) * Lambda.imag();
  Lambda.real() = this->NBSqScratch_->real() * (*this->aointegrals_->ortho2_);
  Lambda.imag() = this->NBSqScratch_->imag() * (*this->aointegrals_->ortho2_);

  (*this->fockA_) += Lambda;
  (*this->fockB_) -= Lambda;

  (*this->fockScalar_) = (*this->fockA_) + (*this->fockB_);
  (*this->fockMz_)     = (*this->fockA_) - (*this->fockB_);
};
*/


template<>
void SingleSlater<dcomplex>::doImagTimeProp(double dt){
  /* Propagate MO coefficients in imaginary time as an alternative to SCF
   * C(new) = exp(-dt * F(current)) * C(current) 
   */
  this->formFock(); // Need orthonormal Fock to propagate
  this->orthoFock3();

  ComplexMatrix propagator = ( -dt * (*this->fockOrthoScalar_) ).exp();
  ComplexMatrix newMOs     = propagator * (*this->moA_);
  *this->moA_              = newMOs; // New MO coefficients are not orthogonal
  propagator               = (*this->moA_).householderQr().householderQ();
  *this->moA_              = propagator;
  if(!this->isClosedShell && this->nTCS_ == 1){
    propagator  = ( -dt * (*this->fockOrthoMz_) ).exp();
    newMOs      = propagator * (*this->moB_);
    *this->moB_ = newMOs; // New MO coefficients are not orthogonal
    propagator  = (*this->moB_).householderQr().householderQ();
    *this->moB_ = propagator;
  }
};




} // namespace ChronusQ
