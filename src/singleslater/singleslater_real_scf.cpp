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
void SingleSlater<double>::complexMem(){;};

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
/*
  RealMap X(this->XMem_,this->nTCS_*this->nBasis_,this->nTCS_*this->nBasis_);
  X = (*this->aointegrals_->overlap_).pow(-0.5); // Make this more efficient... FIXME

  if(this->Ref_ == CUHF){
    RealMap Xp(this->XpMem_,this->nBasis_,this->nBasis_);
    Xp = (*this->aointegrals_->overlap_).pow(0.5); // Make this more efficient... FIXME
  }
*/

  char JOBZ = 'V';
  char UPLO = 'L';
  int INFO;
  auto NTCSxNBASIS = this->nTCS_*this->nBasis_;

  RealVecMap E(this->SEVlMem_,NTCSxNBASIS);
  RealMap    X(this->XMem_   ,NTCSxNBASIS,NTCSxNBASIS);
  RealMap    V(this->SEVcMem_,NTCSxNBASIS,NTCSxNBASIS);
  RealMap    S(this->SCpyMem_,NTCSxNBASIS,NTCSxNBASIS);

  E.setZero();
  V.setZero();
  S.setZero();

  std::memcpy(this->SEVcMem_,this->aointegrals_->overlap_->data(),
    NTCSxNBASIS*NTCSxNBASIS*sizeof(double));

  dsyev_(&JOBZ,&UPLO,&NTCSxNBASIS,this->SEVcMem_,&NTCSxNBASIS,this->SEVlMem_,
    this->WORK_,&this->LWORK_,&INFO);
  
//V.transposeInPlace(); // b/c Row Major...
  std::memcpy(this->SCpyMem_,this->SEVcMem_,NTCSxNBASIS * NTCSxNBASIS *
    sizeof(double));

  for(auto i = 0; i < NTCSxNBASIS; i++)
    S.col(i) /= std::sqrt(this->SEVlMem_[i]);

  X = S * V.transpose();

  if(this->Ref_ == CUHF){
    RealMap    Xp(this->XpMem_    ,NTCSxNBASIS,NTCSxNBASIS);

    for(auto i = 0; i < NTCSxNBASIS; i++)
      S.col(i) *= this->SEVlMem_[i];
 
    Xp = S * V.transpose();
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
//P.transposeInPlace();

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
//FpAlpha.transposeInPlace(); // bc row major
  (*this->moA_) = X * FpAlpha;

  if(!this->isClosedShell && this->Ref_ != TCS){
    FpBeta = X.transpose() * (*this->fockB_) * X;
    dsyev_(&JOBZ,&UPLO,&this->nBasis_,this->FpBetaMem_,&this->nBasis_,this->epsB_->data(),
           this->WORK_,&this->LWORK_,&INFO);
    if(INFO != 0) CErr("DSYEV Failed Fock Beta",this->fileio_->out);
//  FpBeta.transposeInPlace(); // bc row major
    (*this->moB_) = X * FpBeta;
  }

  

}

template<>
void SingleSlater<double>::evalConver(int iter){
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

//if(this->isClosedShell)    this->printDensityInfo(PAlphaRMS,EDelta);
//else if(this->Ref_ != TCS) this->printDensityInfo(PAlphaRMS,PBetaRMS,EDelta);
  if(this->printLevel_ > 0) this->printSCFIter(iter,EDelta,PAlphaRMS,PBetaRMS);
  this->isConverged = (PAlphaRMS < this->denTol_) && (std::pow(EDelta,2) < this->eneTol_);
  if(!this->isClosedShell)
    this->isConverged = this->isConverged && (PBetaRMS < this->denTol_);
  if(this->isPrimary) this->writeSCFFiles();
}

template<>
void SingleSlater<double>::mixOrbitalsSCF(){
  if(this->Ref_ == TCS){
  this->fileio_->out << "** Mixing Alpha-Beta Orbitals for 2C Guess **" << endl;
  //CErr();
  auto nO = this->nAE_ + this->nBE_;
  VectorXd HOMOA,LUMOB;
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

  if(indxHOMOA == -1 || indxLUMOB == -1)
  //  CErr("TCS orbital swap failed to find suitable Alpha-Beta pair",this->fileio_->out);
    return;
  
//CErr();
  HOMOA = this->moA_->col(indxHOMOA) ;
  LUMOB = this->moA_->col(indxLUMOB) ;
//cout << HOMOA << endl << endl;
//cout << LUMOB << endl << endl;
//prettyPrint(cout,*this->moA_,"MO");
  this->moA_->col(indxHOMOA) = std::sqrt(0.5) * (HOMOA + LUMOB);
  this->moA_->col(indxLUMOB) = std::sqrt(0.5) * (HOMOA - LUMOB);
/*
    Eigen::VectorXd HOMO = this->moA_->col(this->nAE_+this->nBE_-1);
    Eigen::VectorXd LUMO = this->moA_->col(this->nTCS_*this->nBasis_-1);
   cout << endl << endl <<  this->moA_->col(this->nAE_+this->nBE_-1) << endl; 
   cout << endl << endl <<  this->moA_->col(this->nTCS_*this->nBasis_-1) << endl;
    this->moA_->col(this->nAE_+this->nBE_-1) = std::sqrt(0.5) * (HOMO + LUMO);
//  this->moA_->col(this->nAE_+this->nBE_) =   std::sqrt(0.5) * (HOMO - LUMO);
    this->moA_->col(this->nTCS_*this->nBasis_-1) = std::sqrt(0.5) * (HOMO - LUMO);

   cout << endl << endl <<  this->moA_->col(this->nAE_+this->nBE_-1) << endl; 
   cout << endl << endl <<  this->moA_->col(this->nTCS_*this->nBasis_-1) << endl;
*/
  }
}

} // namespace ChronusQ
