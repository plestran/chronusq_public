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
using ChronusQ::Controls;
using ChronusQ::FileIO;
using ChronusQ::SingleSlater;

//----------------------------------------//
// do the SCF                             //
// Sajan                                  //
//----------------------------------------//
namespace ChronusQ {
template<>
void SingleSlater<double>::printDensityInfo(double PAlphaRMS,double EDelta){
  this->fileio_->out<<"\nSCF Information:"<<endl;
  this->fileio_->out<<std::right<<std::setw(30)<<"    Delta-E = "<<std::setw(15)<<std::scientific<<EDelta<<std::setw(5)<<" Eh "<<endl;
  this->fileio_->out<<std::right<<std::setw(30)<<"RMS Density = "<<std::setw(15)<<std::scientific<<PAlphaRMS<<endl;
};
template<>
void SingleSlater<double>::printDensityInfo(double PAlphaRMS, double PBetaRMS, double EDelta){
  this->fileio_->out<<"\nSCF Information:"<<endl;
  this->fileio_->out<<std::right<<std::setw(30)<<"    Delta-E = "<<std::setw(15)<<std::scientific<<EDelta<<std::setw(5)<<" Eh "<<endl;
  this->fileio_->out<<std::right<<std::setw(30)<<"RMS Alpha Density = "<<std::setw(15)<<std::scientific<<PAlphaRMS<<endl;
  this->fileio_->out<<std::right<<std::setw(30)<<"RMS Beta Density = "<<std::setw(15)<<std::scientific<<PBetaRMS<<endl;
};

template<>
void SingleSlater<double>::formX(){
  RealMap X(this->XMem_,this->nTCS_*this->nBasis_,this->nTCS_*this->nBasis_);
/*
  if(this->Ref_ != TCS)
    X = (*this->aointegrals_->overlap_).pow(-0.5); // Make this more efficient... FIXME
  else {
    RealMap GenOverlap(this->SMem_,this->nTCS_*this->nBasis_,this->nTCS_*this->nBasis_);
//  RealMatrix OTmp(this->nTCS_*this->nBasis_,this->nTCS_*this->nBasis_);
    for(auto I = 0, i = 0; i < this->nBasis_; i++, I +=2)
    for(auto J = 0, j = 0; j < this->nBasis_; j++, J +=2){
    //OTmp(I,J)    = (*this->aointegrals_->overlap_)(i,j);
    //OTmp(I+1,J+1) = (*this->aointegrals_->overlap_)(i,j);
      GenOverlap(I,J)    = (*this->aointegrals_->overlap_)(i,j);
      GenOverlap(I+1,J+1) = (*this->aointegrals_->overlap_)(i,j);
    }
  //X = OTmp.pow(-0.5);
    X = GenOverlap.pow(-0.5);
  }
*/
  X = (*this->aointegrals_->overlap_).pow(-0.5); // Make this more efficient... FIXME

  if(this->Ref_ == CUHF){
    RealMap Xp(this->XpMem_,this->nBasis_,this->nBasis_);
    Xp = (*this->aointegrals_->overlap_).pow(0.5); // Make this more efficient... FIXME
  }
}

template<>
void SingleSlater<double>::formNO(){
  int INFO;
  char JOBZ = 'V';
  char UPLO = 'L';

  RealMap P(this->PNOMem_,this->nBasis_,this->nBasis_);
  RealMap Xp(this->XpMem_,this->nBasis_,this->nBasis_);

  P = 0.5 * Xp * (*this->densityA_) * Xp;
  if(!this->isClosedShell)
    P += 0.5 * Xp * (*this->densityB_) * Xp;

  dsyev_(&JOBZ,&UPLO,&this->nBasis_,this->PNOMem_,&this->nBasis_,this->occNumMem_,
         this->WORK_,&this->LWORK_,&INFO);
  if(INFO != 0) CErr("DSYEV Failed in FormNO",this->fileio_->out);
  P.transposeInPlace();

  // Swap Ordering
  for(auto i = 0; i < this->nBasis_/2; i++) P.col(i).swap(P.col(this->nBasis_ - i- 1));

}

template<>
void SingleSlater<double>::diagFock(){
  int INFO;
  char JOBZ = 'V';
  char UPLO = 'U';
  auto NTCSxNBASIS = this->nTCS_*this->nBasis_;

  RealMap X(this->XMem_,NTCSxNBASIS,NTCSxNBASIS);
  RealMap POldAlpha(this->POldAlphaMem_,NTCSxNBASIS,NTCSxNBASIS);
  RealMap FpAlpha(this->FpAlphaMem_,NTCSxNBASIS,NTCSxNBASIS);
  RealMap POldBeta(this->POldBetaMem_,0,0);
  RealMap FpBeta(this->FpBetaMem_,0,0);
  if(!this->isClosedShell && this->Ref_ != TCS){
    new (&POldBeta)  RealMap(this->POldBetaMem_, NTCSxNBASIS,NTCSxNBASIS);
    new (&FpBeta)    RealMap(this->FpBetaMem_,NTCSxNBASIS,NTCSxNBASIS);
  }


  if(this->Ref_ == CUHF){
    RealMap P(this->PNOMem_,this->nBasis_,this->nBasis_);
    RealMap Xp(this->XpMem_,this->nBasis_,this->nBasis_);
    RealMap DelF(this->delFMem_,this->nBasis_,this->nBasis_);
    RealMap Lambda(this->lambdaMem_,this->nBasis_,this->nBasis_);

    int activeSpace  = this->molecule_->multip() - 1;
    int coreSpace    = (this->molecule_->nTotalE() - activeSpace) / 2;
    int virtualSpace = this->nBasis_ - coreSpace - activeSpace;

    DelF = 0.5 * X * (*this->fockA_) * X;
    if(!this->isClosedShell)
      DelF -= 0.5 * X * (*this->fockB_) * X;
 
    DelF = P.transpose() * DelF * P;
 
    Lambda.setZero();
    for(auto i = activeSpace + coreSpace; i < this->nBasis_; i++)
    for(auto j = 0                      ; j < coreSpace    ; j++){
      Lambda(i,j) = -DelF(i,j);
      Lambda(j,i) = -DelF(j,i);
    }
    Lambda = P  * Lambda * P.transpose();
    Lambda = Xp * Lambda * Xp;  
 
    (*this->fockA_) += Lambda;
    if(!this->isClosedShell) (*this->fockB_) -= Lambda;
  }

  POldAlpha = (*this->densityA_);
  if(!this->isClosedShell && this->Ref_ != TCS) POldBeta = (*this->densityB_);

  FpAlpha = X.transpose() * (*this->fockA_) * X;
  dsyev_(&JOBZ,&UPLO,&NTCSxNBASIS,this->FpAlphaMem_,&NTCSxNBASIS,this->epsA_->data(),
         this->WORK_,&this->LWORK_,&INFO);
  if(INFO != 0) CErr("DSYEV Failed Fock Alpha",this->fileio_->out);
  FpAlpha.transposeInPlace(); // bc row major
  (*this->moA_) = X * FpAlpha;

  if(!this->isClosedShell && this->Ref_ != TCS){
    FpBeta = X.transpose() * (*this->fockB_) * X;
    dsyev_(&JOBZ,&UPLO,&this->nBasis_,this->FpBetaMem_,&this->nBasis_,this->epsB_->data(),
           this->WORK_,&this->LWORK_,&INFO);
    if(INFO != 0) CErr("DSYEV Failed Fock Beta",this->fileio_->out);
    FpBeta.transposeInPlace(); // bc row major
    (*this->moB_) = X * FpBeta;
  }

  

}

template<>
void SingleSlater<double>::evalConver(){
  double EOld;
  double EDelta;
  double PAlphaRMS;
  double PBetaRMS;

  RealMap POldAlpha(this->POldAlphaMem_,this->nTCS_*this->nBasis_,this->nTCS_*this->nBasis_);
  RealMap POldBeta(this->POldBetaMem_,0,0);
  if(!this->isClosedShell && this->Ref_ != TCS){
    new (&POldBeta)  RealMap(this->POldBetaMem_, this->nBasis_,this->nBasis_);
  }

  EOld = this->totalEnergy;
  this->computeEnergy();
  EDelta = this->totalEnergy - EOld;

  PAlphaRMS = ((*this->densityA_) - POldAlpha).norm();
  if(!this->isClosedShell && this->Ref_ != TCS) PBetaRMS = ((*this->densityB_) - POldBeta).norm();

  if(this->isClosedShell)    this->printDensityInfo(PAlphaRMS,EDelta);
  else if(this->Ref_ != TCS) this->printDensityInfo(PAlphaRMS,PBetaRMS,EDelta);

  this->isConverged = (PAlphaRMS < this->denTol_) && (std::pow(EDelta,2) < this->eneTol_);
  if(!this->isClosedShell)
    this->isConverged = this->isConverged && (PBetaRMS < this->denTol_);
}

template<>
void SingleSlater<double>::mixOrbitalsSCF(){
  if(this->Ref_ == TCS){
    Eigen::VectorXd HOMO = this->moA_->col(this->nAE_+this->nBE_-1);
    Eigen::VectorXd LUMO = this->moA_->col(this->nTCS_*this->nBasis_-1);
    this->moA_->col(this->nAE_+this->nBE_-1) = std::sqrt(0.5) * (HOMO + LUMO);
//  this->moA_->col(this->nAE_+this->nBE_) =   std::sqrt(0.5) * (HOMO - LUMO);
    this->moA_->col(this->nTCS_*this->nBasis_-1) = std::sqrt(0.5) * (HOMO - LUMO);
  }
}

} // namespace ChronusQ
