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
double EDelta;
double PAlphaRMS;
double PBetaRMS;

namespace ChronusQ {
template<>
void SingleSlater<double>::printDensityinf(){
  this->fileio_->out<<"\nSCF Information:"<<endl;
  this->fileio_->out<<std::right<<std::setw(30)<<"    Delta-E = "<<std::setw(15)<<EDelta<<std::setw(5)<<" Eh "<<endl;
  if(this->RHF_)
    this->fileio_->out<<std::right<<std::setw(30)<<"RMS Density = "<<std::setw(15)<<PAlphaRMS<<endl;
  else {
    this->fileio_->out<<std::right<<std::setw(30)<<"RMS Alpha Density = "<<std::setw(15)<<PAlphaRMS<<endl;
    this->fileio_->out<<std::right<<std::setw(30)<<"RMS Beta Density = "<<std::setw(15)<<PBetaRMS<<endl;
  }
};
template<>
void SingleSlater<double>::SCF(){
  if(!this->aointegrals_->haveAOOneE) this->aointegrals_->computeAOOneE();
  double EOld;
  int maxIter    = 128; 
  int n = this->nBasis_; 
  double Dtol = 1e-10;
  double Etol = 1e-8;
  std::vector<RealMatrix> ErrorAlpha;
  std::vector<RealMatrix> fockAlpha;
  std::vector<RealMatrix> ErrorBeta;
  std::vector<RealMatrix> fockBeta;
  int iter; 
  
  //RealMatrix memory allocation
  int lenX = n*n; // X
  int lenF = n*n; // Fock
  int lenP = n*n; // Density
  int lenB = 49;  // B
  int lenCoeff = 7;   // Coeff
  int lwork = 4*n; // LAPACK Workspace
  int LenScr = lenX + lenF + lenP + lenB + lenCoeff + lwork;
  if(!this->RHF_) LenScr += lenF + lenP;

  double *SCR          = NULL;
  double *XMem         = NULL;
  double *FpAlphaMem   = NULL;
  double *POldAlphaMem = NULL;
  double *FpBetaMem    = NULL;
  double *POldBetaMem  = NULL;
  double *BMem         = NULL;
  double *coef         = NULL;
  double *work         = NULL;

  SCR = new double [LenScr];
  XMem    = SCR;
  if(this->RHF_) {
    FpAlphaMem   = XMem + lenX;
    POldAlphaMem = FpAlphaMem + lenF;
    BMem    = POldAlphaMem + lenP;
  } else {
    FpAlphaMem   = XMem + lenX;
    FpBetaMem    = FpAlphaMem + lenF;
    POldAlphaMem = FpBetaMem  + lenF;
    POldBetaMem  = POldAlphaMem + lenP;
    BMem    = POldBetaMem + lenP;
  }

  RealMap X(XMem,n,n);
  RealMap FpAlpha(FpAlphaMem,n,n);
  RealMap POldAlpha(POldAlphaMem,n,n);
  RealMap FpBeta(FpBetaMem,n,n);
  RealMap POldBeta(POldBetaMem,n,n);
  RealMap B(BMem,7,7);
  

  //lapack variables for DIIS
  coef = BMem + lenB;
  int *iPiv = new int[7];
  int row=7;
  int nrhs=1;
  int info=-1;

  //lapack variables for F'C=CE
  char j='V';
  char u='U';
  work = coef + lenCoeff;
   
  
  X=(*this->aointegrals_->overlap_).pow(-0.5);
  for (iter=0; iter<maxIter; iter++){
    this->fileio_->out << endl << endl << bannerTop <<endl;  
    this->fileio_->out << "SCF iteration:"<< iter+1 <<endl;  
    this->fileio_->out << bannerEnd <<endl;  
    
    POldAlpha = (*this->densityA_);
    FpAlpha   = X.transpose()*(*this->fockA_)*X;
    if(!this->RHF_) {
      POldBeta   = (*this->densityB_);
      FpBeta     = X.transpose()*(*this->fockB_)*X;
    }
    EOld = this->totalEnergy;
    dsyev_(&j,&u,&n, FpAlpha.data(), &n, this->epsA_->data(), work, &lwork, &info);
    FpAlpha.transposeInPlace(); // bc Row Major
    (*this->moA_) = X*FpAlpha;
    if(!this->RHF_) {
      dsyev_(&j,&u,&n, FpBeta.data(), &n, this->epsB_->data(), work, &lwork, &info);
      FpBeta.transposeInPlace(); // bc Row Major
      (*this->moB_) = X*FpBeta;
    }
    this->formDensity();
    this->formFock();
    this->computeEnergy();
    
    
    // DIIS 
    if(iter % 6 ==0){
      ErrorAlpha.push_back((*this->fockA_) * (*this->densityA_) * (*this->aointegrals_->overlap_) -
                      (*this->aointegrals_->overlap_) * (*this->densityA_) * (*this->fockA_));
      fockAlpha.push_back(*this->fockA_);
      if(!this->RHF_){
        ErrorBeta.push_back((*this->fockB_) * (*this->densityB_) * (*this->aointegrals_->overlap_) -
                        (*this->aointegrals_->overlap_) * (*this->densityB_) * (*this->fockB_));
        fockBeta.push_back(*this->fockB_);
      }
    }
    if(iter % 6 ==1){
      ErrorAlpha.push_back((*this->fockA_) * (*this->densityA_) * (*this->aointegrals_->overlap_) -
                      (*this->aointegrals_->overlap_) * (*this->densityA_) * (*this->fockA_));
      fockAlpha.push_back(*this->fockA_);
      if(!this->RHF_){
        ErrorBeta.push_back((*this->fockB_) * (*this->densityB_) * (*this->aointegrals_->overlap_) -
                        (*this->aointegrals_->overlap_) * (*this->densityB_) * (*this->fockB_));
        fockBeta.push_back(*this->fockB_);
      }
    }
    if(iter % 6 ==2){
      ErrorAlpha.push_back((*this->fockA_) * (*this->densityA_) * (*this->aointegrals_->overlap_) - 
                      (*this->aointegrals_->overlap_) * (*this->densityA_) * (*this->fockA_));
      fockAlpha.push_back(*this->fockA_);
      if(!this->RHF_){
        ErrorBeta.push_back((*this->fockB_) * (*this->densityB_) * (*this->aointegrals_->overlap_) -
                        (*this->aointegrals_->overlap_) * (*this->densityB_) * (*this->fockB_));
        fockBeta.push_back(*this->fockB_);
      }
    }
    if(iter % 6 ==3){
      ErrorAlpha.push_back((*this->fockA_) * (*this->densityA_) * (*this->aointegrals_->overlap_) - 
                      (*this->aointegrals_->overlap_) * (*this->densityA_) * (*this->fockA_));
      fockAlpha.push_back(*this->fockA_);
      if(!this->RHF_){
        ErrorBeta.push_back((*this->fockB_) * (*this->densityB_) * (*this->aointegrals_->overlap_) -
                        (*this->aointegrals_->overlap_) * (*this->densityB_) * (*this->fockB_));
        fockBeta.push_back(*this->fockB_);
      }
    }
    if(iter % 6 ==4){
      ErrorAlpha.push_back((*this->fockA_) * (*this->densityA_) * (*this->aointegrals_->overlap_) - 
                      (*this->aointegrals_->overlap_) * (*this->densityA_) * (*this->fockA_));
      fockAlpha.push_back(*this->fockA_);
      if(!this->RHF_){
        ErrorBeta.push_back((*this->fockB_) * (*this->densityB_) * (*this->aointegrals_->overlap_) -
                        (*this->aointegrals_->overlap_) * (*this->densityB_) * (*this->fockB_));
        fockBeta.push_back(*this->fockB_);
      }
    }
    if(iter % 6 ==5){
      ErrorAlpha.push_back((*this->fockA_) * (*this->densityA_) * (*this->aointegrals_->overlap_) - 
                      (*this->aointegrals_->overlap_) * (*this->densityA_) * (*this->fockA_));
      fockAlpha.push_back(*this->fockA_);
      if(!this->RHF_){
        ErrorBeta.push_back((*this->fockB_) * (*this->densityB_) * (*this->aointegrals_->overlap_) -
                        (*this->aointegrals_->overlap_) * (*this->densityB_) * (*this->fockB_));
        fockBeta.push_back(*this->fockB_);
      }
    }
    
    
    if(iter % 6==0 && iter!=0){
      
      for (auto j=0;j<ErrorAlpha.size();j++){
        for (auto k=0; k<=j;k++){
          B(j,k)=(ErrorAlpha[j]*(ErrorAlpha[k].transpose())).trace();
          if(!this->RHF_) B(j,k) += (ErrorBeta[j]*(ErrorBeta[k].transpose())).trace();
	  B(k,j)=B(j,k);
        }
      }
      for (auto l=0;l<6;l++){
         B(6,l)=-1.0;
	 B(l,6)=-1.0;
      }
      B(6,6)=0;
      for (auto k =0;k<6;k++){
        coef[k]=0.0;
      }
      coef[6]=-1.0;
      dgesv_(&row,&nrhs,B.data(),&row, iPiv, coef,&row, &info);
      this->fockA_->setZero();
      if(!this->RHF_) this->fockB_->setZero();
      for(auto j = 0; j < 6; j++) {
        *this->fockA_ += coef[j]*fockAlpha[j];
        if(!this->RHF_) *this->fockB_ += coef[j]*fockBeta[j];
      }
      ErrorAlpha.clear();
      fockAlpha.clear();
      if(!this->RHF_){
        ErrorBeta.clear();
        fockBeta.clear();
      }
    }

    PAlphaRMS=((*this->densityA_)-POldAlpha).norm();
    if(!this->RHF_) PBetaRMS = ((*this->densityB_) - POldBeta).norm();
    EDelta= this->totalEnergy-EOld;
    this->printDensityinf();     
    
    if(this->RHF_){
      if(PAlphaRMS < Dtol && pow((this->totalEnergy-EOld),2)<Etol) break;
    } else {
      if(PAlphaRMS < Dtol && PBetaRMS < Dtol && pow((this->totalEnergy-EOld),2)<Etol) break;
    }
  };
  
  //freeing the memory
  delete [] SCR;
  delete [] iPiv;

  this->fileio_->out <<"\n"<<endl; 
  this->fileio_->out << bannerEnd <<endl;
  this->fileio_->out << "\nRequested convergence on RMS density matrix = " <<std::setw(5)<<Dtol <<"  within  " <<maxIter <<"  cycles."<<endl;
  this->fileio_->out << "Requested convergence on             energy = " <<Etol << endl;
  this->fileio_->out << "\nSCF Done: E(RHF) = "<< this->totalEnergy << "  Eh  after  "<< iter+1 << "  cycles" <<endl;
  this->fileio_->out << bannerEnd <<endl;
}; 
} // namespace ChronusQ
