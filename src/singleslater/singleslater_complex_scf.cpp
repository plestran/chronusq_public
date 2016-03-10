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
using ChronusQ::Controls;
using ChronusQ::FileIO;
using ChronusQ::SingleSlater;

//----------------------------------------//
// do the SCF                             //
// Sajan                                  //
//----------------------------------------//
namespace ChronusQ {
template<>
void SingleSlater<dcomplex>::complexMem(){
  this->lenRealScr_ += this->LRWORK_; // Extra space for RWORK for complex
};
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
void SingleSlater<dcomplex>::formX(){
/*
  ComplexMap X(this->XMem_,this->nTCS_*this->nBasis_,this->nTCS_*this->nBasis_);
  X.real() = (*this->aointegrals_->overlap_).pow(-0.5); // Make this more efficient... FIXME

  if(this->Ref_ == CUHF){
    ComplexMap Xp(this->XpMem_,this->nBasis_,this->nBasis_);
    Xp.real() = (*this->aointegrals_->overlap_).pow(0.5); // Make this more efficient... FIXME
  }
*/
  char JOBZ = 'V';
  char UPLO = 'L';
  int INFO;
  auto NTCSxNBASIS = this->nTCS_*this->nBasis_;

  ComplexMap X(this->XMem_   ,NTCSxNBASIS,NTCSxNBASIS);
  RealVecMap E(this->SEVlMem_,NTCSxNBASIS);
  RealMap    V(this->SEVcMem_,NTCSxNBASIS,NTCSxNBASIS);
  RealMap    S(this->SCpyMem_,NTCSxNBASIS,NTCSxNBASIS);

  E.setZero();
  V.setZero();
  S.setZero();

  std::memcpy(this->SEVcMem_,this->aointegrals_->overlap_->data(),
    NTCSxNBASIS*NTCSxNBASIS*sizeof(double));

  dsyev_(&JOBZ,&UPLO,&NTCSxNBASIS,this->SEVcMem_,&NTCSxNBASIS,this->SEVlMem_,
    this->LowdinWORK_,&this->LWORK_,&INFO);
  
//V.transposeInPlace(); // b/c Row Major...
  std::memcpy(this->SCpyMem_,this->SEVcMem_,NTCSxNBASIS * NTCSxNBASIS *
    sizeof(double));

  for(auto i = 0; i < NTCSxNBASIS; i++)
    S.col(i) /= std::sqrt(this->SEVlMem_[i]);

  X.real() = S * V.transpose();

  if(this->Ref_ == CUHF){
    ComplexMap    Xp(this->XpMem_    ,NTCSxNBASIS,NTCSxNBASIS);

    for(auto i = 0; i < NTCSxNBASIS; i++)
      S.col(i) *= this->SEVlMem_[i];
 
    Xp.real() = S * V.transpose();
  }

}

template<>
void SingleSlater<dcomplex>::formNO(){
  int INFO;
  char JOBZ = 'V';
  char UPLO = 'L';

  ComplexMap P(this->PNOMem_,this->nBasis_,this->nBasis_);
  ComplexMap Xp(this->XpMem_,this->nBasis_,this->nBasis_);

  P = 0.5 * Xp * (*this->densityA_) * Xp;
  if(!this->isClosedShell)
    P += 0.5 * Xp * (*this->densityB_) * Xp;

  zheev_(&JOBZ,&UPLO,&this->nBasis_,this->PNOMem_,&this->nBasis_,this->occNumMem_,
         this->WORK_,&this->LWORK_,this->RWORK_,&INFO);
  if(INFO != 0) CErr("ZHEEV Failed in FormNO",this->fileio_->out);
//P.transposeInPlace(); //bc row major

  // Swap Ordering
  for(auto i = 0; i < this->nBasis_/2; i++) P.col(i).swap(P.col(this->nBasis_ - i- 1));

}

template<>
void SingleSlater<dcomplex>::diagFock(){
  int INFO;
  char JOBZ = 'V';
  char UPLO = 'U';
  auto NTCSxNBASIS = this->nTCS_*this->nBasis_;

  ComplexMap X(this->XMem_,NTCSxNBASIS,NTCSxNBASIS);
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
    ComplexMap Xp(this->XpMem_,this->nBasis_,this->nBasis_);
    ComplexMap DelF(this->delFMem_,this->nBasis_,this->nBasis_);
    ComplexMap Lambda(this->lambdaMem_,this->nBasis_,this->nBasis_);

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
  zheev_(&JOBZ,&UPLO,&NTCSxNBASIS,this->FpAlphaMem_,&NTCSxNBASIS,this->epsA_->data(),
         this->WORK_,&this->LWORK_,this->RWORK_,&INFO);
  if(INFO != 0) CErr("ZHEEV Failed Fock Alpha",this->fileio_->out);
//FpAlpha.transposeInPlace(); // bc row major
  (*this->moA_) = X * FpAlpha;

  if(!this->isClosedShell && this->Ref_ != TCS){
    FpBeta = X.transpose() * (*this->fockB_) * X;
    zheev_(&JOBZ,&UPLO,&this->nBasis_,this->FpBetaMem_,&this->nBasis_,this->epsB_->data(),
           this->WORK_,&this->LWORK_,this->RWORK_,&INFO);
    if(INFO != 0) CErr("ZHEEV Failed Fock Beta",this->fileio_->out);
//  FpBeta.transposeInPlace(); // bc row major
    (*this->moB_) = X * FpBeta;
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
 
    PAlphaRMS = ((*this->densityA_).cwiseAbs() - POldAlpha.cwiseAbs()).norm();
    if(!this->isClosedShell && this->Ref_ != TCS) 
      PBetaRMS = ((*this->densityB_).cwiseAbs() - POldBeta.cwiseAbs()).norm();
 
    if(this->printLevel_ > 0) 
      this->printSCFIter(iter,EDelta,PAlphaRMS,PBetaRMS);
 
    this->isConverged = (PAlphaRMS < this->denTol_) && 
                        (std::pow(EDelta,2) < this->eneTol_);

    if(!this->isClosedShell)
      this->isConverged = this->isConverged && (PBetaRMS < this->denTol_);

    if(this->isPrimary) this->writeSCFFiles();
  }
#ifdef CQ_ENABLE_MPI
  MPI_Bcast(&this->isConverged,1,MPI_LOGICAL,0,MPI_COMM_WORLD);
#endif
}

template<>
void SingleSlater<dcomplex>::mixOrbitalsSCF(){
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
  
//CErr();
//indxHOMOA = 2;
//indxLUMOB = 5;
  HOMOA = this->moA_->col(indxHOMOA) ;
  LUMOB = this->moA_->col(indxLUMOB) ;
//cout << HOMOA << endl << endl;
//cout << LUMOB << endl << endl;
//prettyPrintComplex(cout,*this->moA_,"MO");
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

} // namespace ChronusQ
