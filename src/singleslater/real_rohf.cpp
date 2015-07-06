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
extern "C" void dgeqrf_(int *numR, int *numC, double *A, int *lda, double *tau, double *work, int *lwork, int *info);  
//----------------------------------------//
// do the SCF                             //
// Sajan                                  //
//----------------------------------------//

namespace ChronusQ {
/*
template<>
void SingleSlater<double>::printDensityInfo(double PAlphaRMS,double EDelta){
  this->fileio_->out<<"\nSCF Information:"<<endl;
  this->fileio_->out<<std::right<<std::setw(30)<<"    Delta-E = "<<std::setw(15)<<EDelta<<std::setw(5)<<" Eh "<<endl;
  this->fileio_->out<<std::right<<std::setw(30)<<"RMS Density = "<<std::setw(15)<<PAlphaRMS<<endl;
};
template<>
void SingleSlater<double>::printDensityInfo(double PAlphaRMS, double PBetaRMS, double EDelta){
  this->fileio_->out<<"\nSCF Information:"<<endl;
  this->fileio_->out<<std::right<<std::setw(30)<<"    Delta-E = "<<std::setw(15)<<EDelta<<std::setw(5)<<" Eh "<<endl;
  this->fileio_->out<<std::right<<std::setw(30)<<"RMS Alpha Density = "<<std::setw(15)<<PAlphaRMS<<endl;
  this->fileio_->out<<std::right<<std::setw(30)<<"RMS Beta Density = "<<std::setw(15)<<PBetaRMS<<endl;
};
*/
template<>
void SingleSlater<double>::ROHF(){
  if(!this->aointegrals_->haveAOOneE) this->aointegrals_->computeAOOneE();
  int maxIter = 256;
  int n       = this->nBasis_;
  double Dtol = 1e-10;
  double Etol = 1e-08;
  double Eold = this->totalEnergy;
  int iter;
  
  int lwork   = 4*n;
  int lenP    = n*n;
  int lenDelF = n*n;
  int lenX    = n*n;
  int lenFA   = n*n;
  int lenFB   = n*n;
  int lenOccN = n;
  int lenWork = lwork;
  int lenPsort  = n * n;
  int lenXp   = n * n;
  int lenPoldA = n * n;
  int lenPoldB = n * n;
  int lenLambda= n * n;
  int lenCoeff = 7;
  
  int Scr = lenP + lenDelF + lenX + lenFA + lenFB + lenOccN + lenWork + lenPsort + lenXp + lenPoldA + lenPoldB + lenLambda;
  double *SCR, *Pmem, *Delfmem, *Xmem, *XpM, *FAmem, *FBmem, *occNum, *work, *Psortmem, *PoldAmem, *PoldBmem, *lambdAmem;
  double *ErrorAlphaMem  = NULL;  
  double *ErrorBetaMem   = NULL;
  double *FADIIS         = NULL;
  double *FBDIIS         = NULL;
  
  SCR = new double[Scr];
  Pmem = SCR;
  Delfmem = Pmem    + lenP;
  Xmem    = Delfmem + lenDelF;
  FAmem   = Xmem    + lenX;
  FBmem   = FAmem   + lenFA;  
  RealMap P(Pmem,n,n);
  RealMap DelF(Delfmem, n,n);
  RealMap X(Xmem, n,n);
  RealMap FA(FAmem, n,n);
  RealMap FB(FBmem, n,n);

  //Diagonalize P matrix variables
  char j = 'V';
  char u = 'U';
  int info  = -1;
  occNum      = FBmem  + lenFB;
  work        = occNum + lenOccN;
  Psortmem    = work   + lenWork;
  XpM         = Psortmem + lenPsort;
  PoldAmem    = XpM    + lenXp;
  PoldBmem    = PoldAmem + lenPoldA;
  lambdAmem   = PoldBmem + lenPoldB;
  /*
  ErrorAlphaMem = lambdAmem + lenLambda;
  ErrorBetaMem  = ErrorAlphaMem + lenFA*(lenCoeff-1);
  FADIIS = ErrorBetaMem + lenFA*(lenCoeff-1);
  FBDIIS = FADIIS + lenFA*(lenCoeff-1);
  */

  RealMap Xp(XpM,n,n);
  RealMap Psort(Psortmem,n,n);
  RealMap PoldAlpha(PoldAmem,n,n);
  RealMap PoldBeta(PoldBmem,n,n);
  RealMap lambda(lambdAmem,n,n);

  X= (*this->aointegrals_->overlap_).pow(-0.5);
  Xp= (*this->aointegrals_->overlap_).pow(0.5);
  int numE     = this->molecule_->nTotalE();
  int actSpace = this->molecule_->spin() - 1;
  int coreSpace= (numE - actSpace)/2;
  int virSpace = n - coreSpace - actSpace;
  double PAlphaRMS, PBetaRMS, EDelta;

  this->fileio_->out << endl << "Beginnning ROHF Calculation .........." << endl << endl;
  this->moA_->setZero();
  if (!this->RHF_){
    this->moB_->setZero(); 
  }
  
  this->formDensity();
  this->formFock();
  for (iter = 0; iter< maxIter;iter++){
    this->fileio_->out << endl << endl << bannerTop <<endl;  
    this->fileio_->out << "SCF iteration:"<< iter+1 <<endl;  
    this->fileio_->out << bannerEnd <<endl;  
    if (!this->RHF_){
      P = (*this->densityA_ + *this->densityB_)/2 ;
    }
    else{
      P = *this->densityA_/2;
    }
    P = Xp * P * Xp;
    dsyev_(&j, &u, &n, P.data(), &n, occNum, work, &lwork, &info); //diagonalizing P
    P.transposeInPlace();
    
    for (auto i=0;i<n;i++){
      Psort.col(i)=P.col(n-i-1);
    } 
    P = Psort;
    
    if (!this->RHF_){
      DelF = (*this->fockA_ - *this->fockB_) / 2; 
    }
    else{
      DelF = (*this->fockA_)/2;
    }
    DelF = X * DelF * X;
    DelF = P.transpose() * DelF * P;
    lambda.setZero();
    for (auto i=coreSpace+actSpace ; i < n ; i++){
      for (auto j = 0; j<coreSpace; j++){
	lambda(i,j)= -1.0 * DelF(i,j);
        lambda(j,i)= -1.0 * DelF(j,i);
      };
    };
    lambda = P * lambda * P.transpose();
    lambda = Xp *  lambda * Xp;
   
    *this->fockA_ = *this->fockA_ + lambda;
    if (!this->RHF_){
      *this->fockB_ = *this->fockB_ - lambda;
      FB = X.transpose() * (*this->fockB_) * X;
      dsyev_(&j,&u,&n, FB.data(), &n, occNum, work, &lwork, &info);
      FB.transposeInPlace();
      (*this->moB_)= X * FB;
    }
    FA = X.transpose() * (*this->fockA_) * X;
    dsyev_(&j,&u,&n, FA.data(), &n, occNum, work, &lwork, &info);
    FA.transposeInPlace();
    (*this->moA_)= X * FA;
   
    this->formDensity();
    this->formFock();
    this->computeEnergy();
    
    //DIIS
   /*
    RealMap ErrA(ErrorAlphaMem + (iter%(lenCoeff-1))*lenFA,n,n);
    ErrA = (*this->fockA_) * (*this->densityA_) * (*this->aointegrals_->overlap_);
    ErrA -= (*this->aointegrals_->overlap_) * (*this->densityA_) * (*this->fockA_);
    memcpy(FADIIS+(iter%(lenCoeff-1))*lenFA,this->fockA_->data(),lenFA*sizeof(double));
    if(!this->RHF_){
      RealMap ErrB(ErrorBetaMem + (iter%(lenCoeff-1))*lenFA,n,n);
      ErrB = (*this->fockB_) * (*this->densityB_) * (*this->aointegrals_->overlap_);
      ErrB -= (*this->aointegrals_->overlap_) * (*this->densityB_) * (*this->fockB_);
      memcpy(FBDIIS+(iter%(lenCoeff-1))*lenFA,this->fockB_->data(),lenFA*sizeof(double));
    }
    
    if(iter % (lenCoeff-1) == (lenCoeff-2) && iter != 0) 
      this->CDIIS(lenCoeff,ErrorAlphaMem,FADIIS,ErrorBetaMem,FBDIIS);
    */
    PAlphaRMS = ((*this->densityA_)-PoldAlpha).norm();
    if(!this->RHF_) 
      PBetaRMS = ((*this->densityB_)-PoldBeta).norm();
    
    EDelta= this->totalEnergy-Eold;
    if(this->RHF_) this->printDensityInfo(PAlphaRMS,EDelta);     
    else this->printDensityInfo(PAlphaRMS,PBetaRMS,EDelta);     
    
    if(this->RHF_){
      if(PAlphaRMS < Dtol && pow((this->totalEnergy-Eold),2)<Etol) break;
    
    } else {
      if((PAlphaRMS < Dtol || PBetaRMS < Dtol) && pow((this->totalEnergy-Eold),2)<Etol) break;
    }
    Eold = this->totalEnergy;
    PoldAlpha = *this->densityA_;
    if (!this->RHF_)
      PoldBeta = *this->densityB_;
  
  }
  delete[] SCR;
  if(iter >= maxIter)
    this->fileio_->out << "SCF Failed to converge within maximum number of iterations" << endl << endl;
  this->fileio_->out <<"\n"<<endl; 
  this->fileio_->out << bannerEnd <<endl;
  this->fileio_->out << "\nRequested convergence on RMS density matrix = " <<std::setw(5)<<Dtol <<"  within  " <<maxIter <<"  cycles."<<endl;
  this->fileio_->out << "Requested convergence on             energy = " <<Etol << endl;
  if(maxIter > iter){
    this->fileio_->out << "\nSCF Done: E(ROHF) = "<< this->totalEnergy << "  Eh  after  "<< iter+1 << "  cycles" <<endl;
  }
  this->fileio_->out << bannerEnd <<endl;
};
} //namespace ChronusQ
