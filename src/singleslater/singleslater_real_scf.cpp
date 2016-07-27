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
template<>
void SingleSlater<double>::formNO(){
  int INFO;
  char JOBZ = 'V';
  char UPLO = 'L';

  RealMap P(this->PNOMem_,this->nBasis_,this->nBasis_);

  P = 0.5 * (*this->aointegrals_->ortho2_) * (*this->onePDMA_) * 
    (*this->aointegrals_->ortho2_);

  if(!this->isClosedShell)
    P += 0.5 * (*this->aointegrals_->ortho2_) * (*this->onePDMB_) * 
      (*this->aointegrals_->ortho2_);

  int LWORK  = 5*this->nBasis_;
  double *WORK  = this->memManager_->malloc<double>(LWORK);

  dsyev_(&JOBZ,&UPLO,&this->nBasis_,this->PNOMem_,&this->nBasis_,
      this->occNumMem_,WORK,&LWORK,&INFO);

  this->memManager_->free(WORK,LWORK);

  if(INFO != 0) CErr("DSYEV Failed in FormNO",this->fileio_->out);

  // Swap Ordering
  for(auto i = 0; i < this->nBasis_/2; i++) 
    P.col(i).swap(P.col(this->nBasis_ - i - 1));

}

template<>
void SingleSlater<double>::evalConver(int iter){
  double EOld;
  double EDelta;
  double PAlphaRMS;
  double PBetaRMS;
  RealMap POldAlpha (this->POldAlphaMem_,0,0);
  RealMap POldBeta  (this->POldAlphaMem_,0,0);

  if(getRank() == 0) {
    new (&POldAlpha) RealMap(
      this->POldAlphaMem_,this->nTCS_*this->nBasis_,this->nTCS_*this->nBasis_
    );
    if(!this->isClosedShell && this->nTCS_ == 1){
      new (&POldBeta) RealMap(this->POldBetaMem_, this->nBasis_,this->nBasis_);
    }
 
    EOld = this->totalEnergy;
  }
  this->computeEnergy();
  if(getRank() == 0) {
    EDelta = this->totalEnergy - EOld;
 
    PAlphaRMS = ((*this->onePDMA_) - POldAlpha).norm();
    if(!this->isClosedShell && this->nTCS_ == 1) 
      PBetaRMS = ((*this->onePDMB_) - POldBeta).norm();
 
    if(this->printLevel_ > 0) 
      this->printSCFIter(iter,EDelta,PAlphaRMS,PBetaRMS);
    this->isConverged = (PAlphaRMS < this->denTol_) && 
                        (std::abs(EDelta) < this->eneTol_);
    if(!this->isClosedShell)
      this->isConverged = this->isConverged && (PBetaRMS < this->denTol_);

    this->isConverged = this->isConverged || 
      (std::abs(EDelta) < this->eneTol_*5e-2);

    if(this->isPrimary) this->writeSCFFiles();
  }
#ifdef CQ_ENABLE_MPI
  MPI_Bcast(&this->isConverged,1,MPI_LOGICAL,0,MPI_COMM_WORLD);
#endif
}

template<>
void SingleSlater<double>::mixOrbitalsSCF(){
  if(this->nTCS_ != 2) return;

  this->fileio_->out << 
    "** Mixing Alpha-Beta Orbitals for 2C Guess **" << endl;

  auto nO = this->nAE_ + this->nBE_;
  int indxHOMOA = -1, indxLUMOB = -1;

  auto nOrb = this->nBasis_;
  double maxPercentNonZeroAlpha = 0;

  for(auto i = nO-1; i >= 0; i--){
    auto nNonZeroAlpha = 0;
    for(auto j = 0; j < this->nTCS_*this->nBasis_; j+=2){
      auto aComp = (*this->moA_)(j,i);
      auto bComp = (*this->moA_)(j+1,i);
      if(std::abs(aComp) > 1e-10 && std::abs(bComp) < 1e-10) nNonZeroAlpha++;
    }
    double percentNonZeroAlpha = (double)nNonZeroAlpha/(double)nOrb;
    if(percentNonZeroAlpha > maxPercentNonZeroAlpha){
      maxPercentNonZeroAlpha = percentNonZeroAlpha;
      indxHOMOA = i;
    }
  }

  double maxPercentNonZeroBeta = 0;
  for(auto i = nO; i < this->nTCS_*this->nBasis_; i++){
    auto nNonZeroBeta = 0;
    for(auto j = 1; j < this->nTCS_*this->nBasis_; j+=2){
      auto aComp = (*this->moA_)(j-1,i);
      auto bComp = (*this->moA_)(j,i);
      if(std::abs(bComp) > 1e-6 && std::abs(aComp) < 1e-6) nNonZeroBeta++;
    }
    double percentNonZeroBeta = (double)nNonZeroBeta/(double)nOrb;
    if(percentNonZeroBeta > maxPercentNonZeroBeta){
      maxPercentNonZeroBeta = percentNonZeroBeta;
      indxLUMOB = i;
    }
  }

  if(indxHOMOA == -1 || indxLUMOB == -1) return;
  
  RealVecMap HOMOA(this->memManager_->malloc<double>(
        this->nTCS_*this->nBasis_),this->nTCS_*this->nBasis_);
  RealVecMap LUMOB(this->memManager_->malloc<double>(
        this->nTCS_*this->nBasis_),this->nTCS_*this->nBasis_);

  HOMOA = this->moA_->col(indxHOMOA) ;
  LUMOB = this->moA_->col(indxLUMOB) ;
  this->moA_->col(indxHOMOA) = std::sqrt(0.5) * (HOMOA + LUMOB);
  this->moA_->col(indxLUMOB) = std::sqrt(0.5) * (HOMOA - LUMOB);

  this->memManager_->free(HOMOA.data(),this->nTCS_*this->nBasis_);
  this->memManager_->free(LUMOB.data(),this->nTCS_*this->nBasis_);
}

template<>
void SingleSlater<double>::diagFock2(){
  int INFO;
  char JOBZ = 'V';
  char UPLO = 'U';
  auto NTCSxNBASIS = this->nTCS_*this->nBasis_;

  int LWORK  = 5*NTCSxNBASIS;
  double *WORK  = this->memManager_->malloc<double>(LWORK);

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
  this->memManager_->free(WORK,LWORK);
};

template<>
void SingleSlater<double>::orthoFock(){
  if(this->nTCS_ == 1 && this->isClosedShell){
    // F(A)' = X^\dagger * F(A) * X
    this->NBSqScratch_->noalias() = 
      this->aointegrals_->ortho1_->transpose() * (*this->fockA_);
    this->fockOrthoA_->noalias() = 
      (*this->NBSqScratch_) * (*this->aointegrals_->ortho1_);

  } else {
    // F(Scalar)' = X^\dagger * F(Scalar) * X
    (*this->NBSqScratch_) = 
      this->aointegrals_->ortho1_->transpose() * (*this->fockScalar_);
    (*this->fockOrthoScalar_) = 
      (*this->NBSqScratch_) * (*this->aointegrals_->ortho1_);

    // F(Mz)' = X^\dagger * F(Mz) * X
    (*this->NBSqScratch_) = 
      this->aointegrals_->ortho1_->transpose() * (*this->fockMz_);
    (*this->fockOrthoMz_) = 
      (*this->NBSqScratch_) * (*this->aointegrals_->ortho1_);

    std::vector<std::reference_wrapper<TMap>> toGather;
    toGather.emplace_back(*this->fockOrthoScalar_);
    toGather.emplace_back(*this->fockOrthoMz_);
    if(this->nTCS_ == 1)
      // {F(Scalar),F(Mz)} -> {F(A), F(B)}
      Quantum<double>::spinGather(*this->fockOrthoA_,*this->fockOrthoB_,toGather);

    else {
      // F(Mx)' = X^\dagger * F(Mx) * X
      (*this->NBSqScratch_) = 
        this->aointegrals_->ortho1_->transpose() * (*this->fockMx_);
      (*this->fockOrthoMx_) = 
        (*this->NBSqScratch_) * (*this->aointegrals_->ortho1_);

      // F(My)' = X^\dagger * F(My) * X
      (*this->NBSqScratch_) = 
        this->aointegrals_->ortho1_->transpose() * (*this->fockMy_);
      (*this->fockOrthoMy_) = 
        (*this->NBSqScratch_) * (*this->aointegrals_->ortho1_);

      toGather.emplace_back(*this->fockOrthoMy_);
      toGather.emplace_back(*this->fockOrthoMx_);

      // {F(Scalar), F(Mz), F(Mx). F(My)} -> F
      Quantum<double>::spinGather(*this->fockOrthoA_,toGather);
    }
  }
};

template<>
void SingleSlater<double>::fockCUHF() {
  RealMap P(this->PNOMem_,this->nBasis_,this->nBasis_);
  RealMap DelF(this->delFMem_,this->nBasis_,this->nBasis_);
  RealMap Lambda(this->lambdaMem_,this->nBasis_,this->nBasis_);

  int activeSpace  = this->molecule_->multip() - 1;
  int coreSpace    = (this->molecule_->nTotalE() - activeSpace) / 2;
  int virtualSpace = this->nBasis_ - coreSpace - activeSpace;

  // DelF = X * (F(A) - F(B)) * X
  (*this->NBSqScratch_) = 0.5 * (*this->aointegrals_->ortho1_) *
    (*this->fockMz_);
  DelF = (*this->NBSqScratch_) * (*this->aointegrals_->ortho1_);

  // DelF = C(NO)^\dagger * DelF * C(NO) (Natural Orbitals)
  (*this->NBSqScratch_) = P.transpose() * DelF;
  DelF = (*this->NBSqScratch_) * P;

  Lambda.setZero();
  for(auto i = activeSpace + coreSpace; i < this->nBasis_; i++)
  for(auto j = 0                      ; j < coreSpace    ; j++){
    Lambda(i,j) = -DelF(i,j);
    Lambda(j,i) = -DelF(j,i);
  }

  (*this->NBSqScratch_) = P * Lambda;
  Lambda = (*this->NBSqScratch_) * P.transpose();

  (*this->NBSqScratch_) = (*this->aointegrals_->ortho2_) * Lambda;
  Lambda = (*this->NBSqScratch_) * (*this->aointegrals_->ortho2_);

  (*this->fockA_) += Lambda;
  (*this->fockB_) -= Lambda;

  (*this->fockScalar_) = (*this->fockA_) + (*this->fockB_);
  (*this->fockMz_)     = (*this->fockA_) - (*this->fockB_);
};


template<>
void SingleSlater<double>::orthoDen(){
  if(this->nTCS_ == 1 && this->isClosedShell) {
    (*this->NBSqScratch_) = 
      (*this->aointegrals_->ortho1_) * (*this->onePDMA_);
    (*this->onePDMOrthoA_) = 
      (*this->NBSqScratch_) * (*this->aointegrals_->ortho1_);
    (*this->onePDMA_) = (*this->onePDMOrthoA_);

  } else {
    std::vector<std::reference_wrapper<RealMap>> scattered;
    scattered.emplace_back(*this->onePDMOrthoScalar_);
    scattered.emplace_back(*this->onePDMOrthoMz_);
    if(this->nTCS_ == 1) {
      Quantum<double>::spinScatter(*this->onePDMA_,*this->onePDMB_,scattered);
    } else {
      scattered.emplace_back(*this->onePDMOrthoMy_);
      scattered.emplace_back(*this->onePDMOrthoMx_);
      Quantum<double>::spinScatter(*this->onePDMA_,scattered);
    }

    (*this->NBSqScratch_) = 
      (*this->aointegrals_->ortho1_) * (*this->onePDMOrthoScalar_);
    (*this->onePDMOrthoScalar_) = 
      (*this->NBSqScratch_) * (*this->aointegrals_->ortho1_);
    (*this->onePDMScalar_) = (*this->onePDMOrthoScalar_);

    (*this->NBSqScratch_) = 
      (*this->aointegrals_->ortho1_) * (*this->onePDMOrthoMz_);
    (*this->onePDMOrthoMz_) = 
      (*this->NBSqScratch_) * (*this->aointegrals_->ortho1_);
    (*this->onePDMMz_) = (*this->onePDMOrthoMz_);

    std::vector<std::reference_wrapper<RealMap>> toGather;
    toGather.emplace_back(*this->onePDMScalar_);
    toGather.emplace_back(*this->onePDMMz_);

    if(this->nTCS_ == 2) {
      (*this->NBSqScratch_) = 
        (*this->aointegrals_->ortho1_) * (*this->onePDMOrthoMx_);
      (*this->onePDMOrthoMx_) = 
        (*this->NBSqScratch_) * (*this->aointegrals_->ortho1_);
      (*this->onePDMMx_) = (*this->onePDMOrthoMx_);

      (*this->NBSqScratch_) = 
        (*this->aointegrals_->ortho1_) * (*this->onePDMOrthoMy_);
      (*this->onePDMOrthoMy_) = 
        (*this->NBSqScratch_) * (*this->aointegrals_->ortho1_);
      (*this->onePDMMy_) = (*this->onePDMOrthoMy_);

      toGather.emplace_back(*this->onePDMMy_);
      toGather.emplace_back(*this->onePDMMx_);
      Quantum<double>::spinGather(*this->onePDMA_,toGather);
    } else
      Quantum<double>::spinGather(*this->onePDMA_,*this->onePDMB_,toGather);
  }
};

} // namespace ChronusQ
