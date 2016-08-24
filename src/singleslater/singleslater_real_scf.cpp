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

namespace ChronusQ {
template<>
void SingleSlater<double>::formNO(){
  int INFO;
  char JOBZ = 'V';
  char UPLO = 'L';

  RealMap P(this->PNOMem_,this->nBasis_,this->nBasis_);

/*
  P = 0.5 * (*this->aointegrals_->ortho2_) * (*this->onePDMA_) * 
    (*this->aointegrals_->ortho2_);

  if(!this->isClosedShell)
    P += 0.5 * (*this->aointegrals_->ortho2_) * (*this->onePDMB_) * 
      (*this->aointegrals_->ortho2_);
*/
  this->aointegrals_->Ortho2Trans(*this->onePDMScalar_,P);
  P *= 0.5;

  int LWORK  = 5*this->nBasis_;
  double *WORK  = this->memManager_->malloc<double>(LWORK);

  dsyev_(&JOBZ,&UPLO,&this->nBasis_,this->PNOMem_,&this->nBasis_,
      this->occNumMem_,WORK,&LWORK,&INFO);
  RealMap OCC(this->occNumMem_,this->nBasis_,1);
  prettyPrint(cout,OCC,"OCC");

  this->memManager_->free(WORK,LWORK);

  if(INFO != 0) CErr("DSYEV Failed in FormNO",this->fileio_->out);

  // Swap Ordering
  for(auto i = 0; i < this->nBasis_/2; i++) 
    P.col(i).swap(P.col(this->nBasis_ - i - 1));

}


template<>
void SingleSlater<double>::mixOrbitalsComplex(){ };

template<>
void SingleSlater<double>::diagFock2(){
  int INFO;
  char JOBZ = 'V';
  char UPLO = 'U';
  auto NTCSxNBASIS = this->nTCS_*this->nBasis_;

  int LWORK  = 5*NTCSxNBASIS;
  double *WORK  = this->memManager_->malloc<double>(LWORK);

/*
  dsyev_(&JOBZ,&UPLO,&NTCSxNBASIS,this->fockOrthoA_->data(),&NTCSxNBASIS,
      this->epsA_->data(),WORK,&LWORK,&INFO);
  if(INFO != 0) CErr("DSYEV Failed Fock Alpha",this->fileio_->out);
  (*this->moA_) = (*this->fockOrthoA_);

  if(this->nTCS_ == 1 && !this->isClosedShell){
    dsyev_(&JOBZ,&UPLO,&NTCSxNBASIS,this->fockOrthoB_->data(),&NTCSxNBASIS,
        this->epsB_->data(),WORK,&LWORK,&INFO);
    if(INFO != 0) CErr("DSYEV Failed Fock Beta",this->fileio_->out);
    (*this->moB_) = (*this->fockOrthoB_);
  }
*/
  // ** WARNING **
  // this assumes that the orthonormal fock has been copied into
  // the MO coefficients storage
  dsyev_(&JOBZ,&UPLO,&NTCSxNBASIS,this->moA_->data(),&NTCSxNBASIS,
      this->epsA_->data(),WORK,&LWORK,&INFO);
  if(INFO != 0) CErr("DSYEV Failed Fock Alpha",this->fileio_->out);
  prettyPrintSmart(cout,*this->moA_,"MOA");

  if(this->nTCS_ == 1 && !this->isClosedShell){
    dsyev_(&JOBZ,&UPLO,&NTCSxNBASIS,this->moB_->data(),&NTCSxNBASIS,
        this->epsB_->data(),WORK,&LWORK,&INFO);
    if(INFO != 0) CErr("DSYEV Failed Fock Beta",this->fileio_->out);
    prettyPrintSmart(cout,*this->moB_,"MOB");
  }
  this->memManager_->free(WORK,LWORK);
};


template<>
void SingleSlater<double>::doImagTimeProp(double dt){
  /* Propagate MO coefficients in imaginary time as an alternative to SCF
   * C(new) = exp(-dt * F(current)) * C(current) 
   */
  this->formFock(); // Need orthonormal Fock to propagate
  this->orthoFock3();

  RealMatrix propagator = ( -dt * (*this->fockOrthoScalar_) ).exp();
  RealMatrix newMOs     = propagator * (*this->moA_);
  *this->moA_           = newMOs; // New MO coefficients are not orthogonal
  propagator            = (*this->moA_).householderQr().householderQ();
  *this->moA_           = propagator;
  if(!this->isClosedShell && this->nTCS_ == 1){
    propagator  = ( -dt * (*this->fockOrthoMz_) ).exp();
    newMOs      = propagator * (*this->moB_);
    *this->moB_ = newMOs; // New MO coefficients are not orthogonal
    propagator  = (*this->moB_).householderQr().householderQ();
    *this->moB_ = propagator;
    }
};



} // namespace ChronusQ
