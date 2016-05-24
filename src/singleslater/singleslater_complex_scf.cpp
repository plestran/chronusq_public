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
void SingleSlater<dcomplex>::printDensityInfo(double PAlphaRMS,double EDelta){
  this->fileio_->out<<"\nSCF Information:"<<endl;
  this->fileio_->out<<std::right<<std::setw(30)<<"    Delta-E = "<<std::setw(15)<<std::scientific<<EDelta<<std::setw(5)<<" Eh "<<endl;
  this->fileio_->out<<std::right<<std::setw(30)<<"RMS Density = "<<std::setw(15)<<std::scientific<<PAlphaRMS<<endl;
};
template<>
void SingleSlater<dcomplex>::printDensityInfo(double PAlphaRMS, double PBetaRMS, double EDelta){
  this->fileio_->out<<"\nSCF Information:"<<endl;
  this->fileio_->out<<std::right<<std::setw(30)<<"    Delta-E = "<<std::setw(15)<<std::scientific<<EDelta<<std::setw(5)<<" Eh "<<endl;
  this->fileio_->out<<std::right<<std::setw(30)<<"RMS Alpha Density = "<<std::setw(15)<<std::scientific<<PAlphaRMS<<endl;
  this->fileio_->out<<std::right<<std::setw(30)<<"RMS Beta Density = "<<std::setw(15)<<std::scientific<<PBetaRMS<<endl;
};


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

  zheev_(&JOBZ,&UPLO,&this->nBasis_,this->PNOMem_,&this->nBasis_,
      this->occNumMem_,this->WORK_,&this->LWORK_,this->RWORK_,&INFO);

  if(INFO != 0) CErr("ZHEEV Failed in FormNO",this->fileio_->out);
//P.transposeInPlace(); //bc row major

  // Swap Ordering
  for(auto i = 0; i < this->nBasis_/2; i++) 
    P.col(i).swap(P.col(this->nBasis_ - i- 1));

}

template<>
void SingleSlater<dcomplex>::diagFock(){
  int INFO;
  char JOBZ = 'V';
  char UPLO = 'U';
  auto NTCSxNBASIS = this->nTCS_*this->nBasis_;

  ComplexMap POldAlpha(this->POldAlphaMem_,NTCSxNBASIS,NTCSxNBASIS);
  ComplexMap FpAlpha(this->FpAlphaMem_,NTCSxNBASIS,NTCSxNBASIS);
  ComplexMap POldBeta(this->POldBetaMem_,0,0);
  ComplexMap FpBeta(this->FpBetaMem_,0,0);
  if(!this->isClosedShell && this->Ref_ != TCS){
    new (&POldBeta)  ComplexMap(this->POldBetaMem_, NTCSxNBASIS,NTCSxNBASIS);
    new (&FpBeta)    ComplexMap(this->FpBetaMem_,NTCSxNBASIS,NTCSxNBASIS);
  }


  if(this->Ref_ == CUHF){
    ComplexMap P(this->PNOMem_,this->nBasis_,this->nBasis_);
    ComplexMap DelF(this->delFMem_,this->nBasis_,this->nBasis_);
    ComplexMap Lambda(this->lambdaMem_,this->nBasis_,this->nBasis_);

    int activeSpace  = this->molecule_->multip() - 1;
    int coreSpace    = (this->molecule_->nTotalE() - activeSpace) / 2;
    int virtualSpace = this->nBasis_ - coreSpace - activeSpace;

    DelF.real() = 0.5 * (*this->aointegrals_->ortho1_) * this->fockA_->real() * 
      (*this->aointegrals_->ortho1_);
    DelF.imag() = 0.5 * (*this->aointegrals_->ortho1_) * this->fockA_->imag() * 
      (*this->aointegrals_->ortho1_);

    if(!this->isClosedShell){
      DelF.real() -= 0.5 * (*this->aointegrals_->ortho1_) * 
        this->fockB_->real() * (*this->aointegrals_->ortho1_);
      DelF.imag() -= 0.5 * (*this->aointegrals_->ortho1_) * 
        this->fockB_->imag() * (*this->aointegrals_->ortho1_);
    }
 
    DelF = P.transpose() * DelF * P;
 
    Lambda.setZero();
    for(auto i = activeSpace + coreSpace; i < this->nBasis_; i++)
    for(auto j = 0                      ; j < coreSpace    ; j++){
      Lambda(i,j) = -DelF(i,j);
      Lambda(j,i) = -DelF(j,i);
    }

    Lambda = P  * Lambda * P.transpose();
    Lambda.real() = (*this->aointegrals_->ortho2_) * Lambda.real() * 
      (*this->aointegrals_->ortho2_);  
    Lambda.imag() = (*this->aointegrals_->ortho2_) * Lambda.imag() * 
      (*this->aointegrals_->ortho2_);  
 
    (*this->fockA_) += Lambda;
    if(!this->isClosedShell) (*this->fockB_) -= Lambda;
  }

  POldAlpha = (*this->onePDMA_);
  if(!this->isClosedShell && this->Ref_ != TCS) POldBeta = (*this->onePDMB_);

  FpAlpha.real() = (*this->aointegrals_->ortho1_).transpose() * 
    this->fockA_->real() * (*this->aointegrals_->ortho1_);
  FpAlpha.imag() = (*this->aointegrals_->ortho1_).transpose() * 
    this->fockA_->imag() * (*this->aointegrals_->ortho1_);

  zheev_(&JOBZ,&UPLO,&NTCSxNBASIS,this->FpAlphaMem_,&NTCSxNBASIS,
      this->epsA_->data(),this->WORK_,&this->LWORK_,this->RWORK_,&INFO);

  if(INFO != 0) CErr("ZHEEV Failed Fock Alpha",this->fileio_->out);
//FpAlpha.transposeInPlace(); // bc row major
  this->moA_->real() = (*this->aointegrals_->ortho1_) * FpAlpha.real();
  this->moA_->imag() = (*this->aointegrals_->ortho1_) * FpAlpha.imag();

  if(!this->isClosedShell && this->Ref_ != TCS){
    FpBeta.real() = (*this->aointegrals_->ortho1_).transpose() * 
      this->fockB_->real() * (*this->aointegrals_->ortho1_);
    FpBeta.imag() = (*this->aointegrals_->ortho1_).transpose() * 
      this->fockB_->imag() * (*this->aointegrals_->ortho1_);

    zheev_(&JOBZ,&UPLO,&this->nBasis_,this->FpBetaMem_,&this->nBasis_,
        this->epsB_->data(),this->WORK_,&this->LWORK_,this->RWORK_,&INFO);

    if(INFO != 0) CErr("ZHEEV Failed Fock Beta",this->fileio_->out);
//  FpBeta.transposeInPlace(); // bc row major
    this->moB_->real() = (*this->aointegrals_->ortho1_) * FpBeta.real();
    this->moB_->imag() = (*this->aointegrals_->ortho1_) * FpBeta.imag();
  }

  

}

template<>
void SingleSlater<dcomplex>::evalConver(int iter){
  double EOld;
  double EDelta;
  double PAlphaRMS;
  double PBetaRMS;
  ComplexMap POldAlpha (this->POldAlphaMem_,0,0);
  ComplexMap POldBeta  (this->POldAlphaMem_,0,0);

  if(getRank() == 0){
    new (&POldAlpha) ComplexMap(
      this->POldAlphaMem_,this->nTCS_*this->nBasis_,this->nTCS_*this->nBasis_
    );
    if(!this->isClosedShell && this->Ref_ != TCS){
      new (&POldBeta) ComplexMap(
        this->POldBetaMem_,this->nBasis_,this->nBasis_
      );
    }
 
    EOld = this->totalEnergy;
  }
  this->computeEnergy();
  if(getRank() == 0){
    EDelta = this->totalEnergy - EOld;
 
    PAlphaRMS = ((*this->onePDMA_).cwiseAbs() - POldAlpha.cwiseAbs()).norm();
    if(!this->isClosedShell && this->Ref_ != TCS) 
      PBetaRMS = ((*this->onePDMB_).cwiseAbs() - POldBeta.cwiseAbs()).norm();
 
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
void SingleSlater<dcomplex>::mixOrbitalsSCF(){
  return;
  auto nO = this->nAE_ + this->nBE_;
  if(this->Ref_ == TCS){
  //CErr();
  this->fileio_->out << "** Mixing Alpha-Beta Orbitals for 2C Guess **" << endl;
  Eigen::VectorXcd HOMOA,LUMOB;
  int indxHOMOA = -1, indxLUMOB = -1;
/*
  for(auto i = nO-1; i >= 0; i--){
    auto aComp = this->moA_->col(i)(0);
    auto bComp = this->moA_->col(i)(1);
    if(std::abs(aComp) > 1e-10 && std::abs(bComp) < 1e-10){
      HOMOA = this->moA_->col(i);
      indxHOMOA = i;
      break;
    }
  }
  for(auto i = nO; i < this->nTCS_*this->nBasis_; i++){
    auto aComp = this->moA_->col(i)(0);
    auto bComp = this->moA_->col(i)(1);
    if(std::abs(bComp) > 1e-10 && std::abs(aComp) < 1e-10){
      LUMOB = this->moA_->col(i);
      indxLUMOB = i;
      break;
    }
  }
*/
  auto nOrb = this->nBasis_;
  double maxPercentNonZeroAlpha = 0;
  for(auto i = nO-1; i >= 0; i--){
    auto nNonZeroAlpha = 0;
    for(auto j = 0; j < this->nTCS_*this->nBasis_; j+=2){
      auto aComp = (*this->moA_)(j,i);
      auto bComp = (*this->moA_)(j+1,i);
      if(std::norm(aComp) > 1e-12 && std::norm(bComp) < 1e-12) nNonZeroAlpha++;
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
      if(std::norm(bComp) > 1e-12 && std::norm(aComp) < 1e-12) nNonZeroBeta++;
    }
    double percentNonZeroBeta = (double)nNonZeroBeta/(double)nOrb;
    if(percentNonZeroBeta > maxPercentNonZeroBeta){
      maxPercentNonZeroBeta = percentNonZeroBeta;
      indxLUMOB = i;
    }
  }

  if(indxHOMOA == -1 || indxLUMOB == -1){
    this->fileio_->out << "TCS orbital swap failed to find suitable Alpha-Beta pair" << endl;
    return;
  }
  
  HOMOA = this->moA_->col(indxHOMOA) ;
  LUMOB = this->moA_->col(indxLUMOB) ;
  this->moA_->col(indxHOMOA) = std::sqrt(0.5) * (HOMOA + LUMOB);
  this->moA_->col(indxLUMOB) = std::sqrt(0.5) * (HOMOA - LUMOB);
  }
  this->fileio_->out << "** Mixing HOMO and LUMO for Complex Guess **" << endl;
  if (this->Ref_==TCS) {
    auto HOMO = this->moA_->col(nO-1);
    auto LUMO = this->moA_->col(nO);
    this->moA_->col(nO-1) = std::sqrt(0.5) * (HOMO + math.ii*LUMO);
    this->moA_->col(nO)   = std::sqrt(0.5) * (HOMO - math.ii*LUMO);
  } else {
    auto HOMO = this->moA_->col(this->nOccA_-1);
    auto LUMO = this->moA_->col(this->nOccA_);
    this->moA_->col(this->nOccA_-1) = std::sqrt(0.5) * (HOMO + math.ii*LUMO);
    this->moA_->col(this->nOccA_)   = std::sqrt(0.5) * (HOMO - math.ii*LUMO);
  }
}

template<>
void SingleSlater<dcomplex>::diagFock2(){
  int INFO;
  char JOBZ = 'V';
  char UPLO = 'U';
  auto NTCSxNBASIS = this->nTCS_*this->nBasis_;

  zheev_(&JOBZ,&UPLO,&NTCSxNBASIS,this->fockOrthoA_->data(),&NTCSxNBASIS,
      this->epsA_->data(),this->WORK_,&this->LWORK_,this->RWORK_,&INFO);
  if(INFO != 0) CErr("ZHEEV Failed Fock Alpha",this->fileio_->out);
  (*this->moA_) = (*this->fockOrthoA_);

  if(this->nTCS_ == 1 && !this->isClosedShell){
    zheev_(&JOBZ,&UPLO,&NTCSxNBASIS,this->fockOrthoB_->data(),&NTCSxNBASIS,
        this->epsB_->data(),this->WORK_,&this->LWORK_,this->RWORK_,&INFO);
    if(INFO != 0) CErr("ZHEEV Failed Fock Beta",this->fileio_->out);
    (*this->moB_) = (*this->fockOrthoB_);
  }
};

template<>
void SingleSlater<dcomplex>::orthoFock(){
  if(this->nTCS_ == 1 && this->isClosedShell){
    // F(A)' = X^\dagger * F(A) * X
    this->NBSqScratch_->real() = 
      this->aointegrals_->ortho1_->transpose() * this->fockA_->real();
    this->NBSqScratch_->imag() = 
      this->aointegrals_->ortho1_->transpose() * this->fockA_->imag();

    this->fockOrthoA_->real() = 
      this->NBSqScratch_->real() * (*this->aointegrals_->ortho1_);
    this->fockOrthoA_->imag() = 
      this->NBSqScratch_->imag() * (*this->aointegrals_->ortho1_);


  } else {
    // F(Scalar)' = X^\dagger * F(Scalar) * X
    this->NBSqScratch_->real() = 
      this->aointegrals_->ortho1_->transpose() * this->fockScalar_->real();
    this->NBSqScratch_->imag() = 
      this->aointegrals_->ortho1_->transpose() * this->fockScalar_->imag();

    this->fockOrthoScalar_->real() = 
      this->NBSqScratch_->real() * (*this->aointegrals_->ortho1_);
    this->fockOrthoScalar_->imag() = 
      this->NBSqScratch_->imag() * (*this->aointegrals_->ortho1_);

    // F(Mz)' = X^\dagger * F(Mz) * X
    this->NBSqScratch_->real() = 
      this->aointegrals_->ortho1_->transpose() * this->fockMz_->real();
    this->NBSqScratch_->imag() = 
      this->aointegrals_->ortho1_->transpose() * this->fockMz_->imag();

    this->fockOrthoMz_->real() = 
      this->NBSqScratch_->real() * (*this->aointegrals_->ortho1_);
    this->fockOrthoMz_->imag() = 
      this->NBSqScratch_->imag() * (*this->aointegrals_->ortho1_);

    std::vector<std::reference_wrapper<TMap>> toGather;
    toGather.emplace_back(*this->fockOrthoScalar_);
    toGather.emplace_back(*this->fockOrthoMz_);
    if(this->nTCS_ == 1)
      // {F(Scalar),F(Mz)} -> {F(A), F(B)}
      Quantum<dcomplex>::spinGather(*this->fockOrthoA_,*this->fockOrthoB_,
          toGather);
    else {
      // F(Mx)' = X^\dagger * F(Mx) * X
      this->NBSqScratch_->real() = 
        this->aointegrals_->ortho1_->transpose() * this->fockMx_->real();
      this->NBSqScratch_->imag() = 
        this->aointegrals_->ortho1_->transpose() * this->fockMx_->imag();
     
      this->fockOrthoMx_->real() = 
        this->NBSqScratch_->real() * (*this->aointegrals_->ortho1_);
      this->fockOrthoMx_->imag() = 
        this->NBSqScratch_->imag() * (*this->aointegrals_->ortho1_);

      // F(My)' = X^\dagger * F(My) * X
      this->NBSqScratch_->real() = 
        this->aointegrals_->ortho1_->transpose() * this->fockMy_->real();
      this->NBSqScratch_->imag() = 
        this->aointegrals_->ortho1_->transpose() * this->fockMy_->imag();
     
      this->fockOrthoMy_->real() = 
        this->NBSqScratch_->real() * (*this->aointegrals_->ortho1_);
      this->fockOrthoMy_->imag() = 
        this->NBSqScratch_->imag() * (*this->aointegrals_->ortho1_);

      toGather.emplace_back(*this->fockOrthoMy_);
      toGather.emplace_back(*this->fockOrthoMx_);

      // {F(Scalar), F(Mz), F(Mx). F(My)} -> F
      Quantum<dcomplex>::spinGather(*this->fockOrthoA_,toGather);
    }
  }
};

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

template<>
void SingleSlater<dcomplex>::orthoDen(){
  if(this->nTCS_ == 1 && this->isClosedShell) {
    this->NBSqScratch_->real() = 
      (*this->aointegrals_->ortho1_) * this->onePDMA_->real();
    this->NBSqScratch_->imag() = 
      (*this->aointegrals_->ortho1_) * this->onePDMA_->imag();

    this->onePDMOrthoA_->real() = 
      this->NBSqScratch_->real() * (*this->aointegrals_->ortho1_);
    this->onePDMOrthoA_->imag() = 
      this->NBSqScratch_->imag() * (*this->aointegrals_->ortho1_);

    (*this->onePDMA_) = (*this->onePDMOrthoA_);

  } else {
    std::vector<std::reference_wrapper<ComplexMap>> scattered;
    scattered.emplace_back(*this->onePDMOrthoScalar_);
    scattered.emplace_back(*this->onePDMOrthoMz_);
    if(this->nTCS_ == 1) {
      Quantum<dcomplex>::spinScatter(*this->onePDMA_,*this->onePDMB_,scattered);
    } else {
      scattered.emplace_back(*this->onePDMOrthoMy_);
      scattered.emplace_back(*this->onePDMOrthoMx_);
      Quantum<dcomplex>::spinScatter(*this->onePDMA_,scattered);
    }

    this->NBSqScratch_->real() = 
      (*this->aointegrals_->ortho1_) * this->onePDMScalar_->real();
    this->NBSqScratch_->imag() = 
      (*this->aointegrals_->ortho1_) * this->onePDMScalar_->imag();

    this->onePDMOrthoScalar_->real() = 
      this->NBSqScratch_->real() * (*this->aointegrals_->ortho1_);
    this->onePDMOrthoScalar_->imag() = 
      this->NBSqScratch_->imag() * (*this->aointegrals_->ortho1_);

    (*this->onePDMScalar_) = (*this->onePDMOrthoScalar_);

    this->NBSqScratch_->real() = 
      (*this->aointegrals_->ortho1_) * this->onePDMMz_->real();
    this->NBSqScratch_->imag() = 
      (*this->aointegrals_->ortho1_) * this->onePDMMz_->imag();

    this->onePDMOrthoMz_->real() = 
      this->NBSqScratch_->real() * (*this->aointegrals_->ortho1_);
    this->onePDMOrthoMz_->imag() = 
      this->NBSqScratch_->imag() * (*this->aointegrals_->ortho1_);

    (*this->onePDMMz_) = (*this->onePDMOrthoMz_);

    std::vector<std::reference_wrapper<ComplexMap>> toGather;
    toGather.emplace_back(*this->onePDMScalar_);
    toGather.emplace_back(*this->onePDMMz_);

    if(this->nTCS_ == 2) {
      this->NBSqScratch_->real() = 
        (*this->aointegrals_->ortho1_) * this->onePDMMx_->real();
      this->NBSqScratch_->imag() = 
        (*this->aointegrals_->ortho1_) * this->onePDMMx_->imag();
     
      this->onePDMOrthoMx_->real() = 
        this->NBSqScratch_->real() * (*this->aointegrals_->ortho1_);
      this->onePDMOrthoMx_->imag() = 
        this->NBSqScratch_->imag() * (*this->aointegrals_->ortho1_);

      (*this->onePDMMx_) = (*this->onePDMOrthoMx_);

      this->NBSqScratch_->real() = 
        (*this->aointegrals_->ortho1_) * this->onePDMMy_->real();
      this->NBSqScratch_->imag() = 
        (*this->aointegrals_->ortho1_) * this->onePDMMy_->imag();
     
      this->onePDMOrthoMy_->real() = 
        this->NBSqScratch_->real() * (*this->aointegrals_->ortho1_);
      this->onePDMOrthoMy_->imag() = 
        this->NBSqScratch_->imag() * (*this->aointegrals_->ortho1_);

      (*this->onePDMMy_) = (*this->onePDMOrthoMy_);

      toGather.emplace_back(*this->onePDMMy_);
      toGather.emplace_back(*this->onePDMMx_);
      Quantum<dcomplex>::spinGather(*this->onePDMA_,toGather);
    } else
      Quantum<dcomplex>::spinGather(*this->onePDMA_,*this->onePDMB_,toGather);
  }
};

} // namespace ChronusQ
