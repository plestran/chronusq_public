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
template<>
void SingleSlater<double>::SCF(){
  if(!this->aointegrals_->haveAOOneE) this->aointegrals_->computeAOOneE();
  double EOld;
  double EDelta;
  double PAlphaRMS;
  double PBetaRMS;

  int maxIter    = 128; 
  int n = this->nBasis_; 
  double Dtol = 1e-10;
  double Etol = 1e-8;
  int iter; 
  
  //RealMatrix memory allocation
  int lenX = n*n; // X
  int lenF = n*n; // Fock
  int lenP = n*n; // Density
  int lenB = 49;  // B
  int lenCoeff = 7;   // Coeff
  int lwork = 4*n; // LAPACK Workspace
  int LenScr = lenX + lenF + lenP + lenB + lenCoeff + lwork + 2*lenF*(lenCoeff-1);
  if(!this->RHF_) LenScr += lenF + lenP + 2*lenF*(lenCoeff-1);

  double *SCR            = NULL;
  double *XMem           = NULL;
  double *FpAlphaMem     = NULL;
  double *POldAlphaMem   = NULL;
  double *FpBetaMem      = NULL;
  double *POldBetaMem    = NULL;
  double *BMem           = NULL;
  double *coef           = NULL;
  double *ErrorAlphaMem  = NULL;  
  double *ErrorBetaMem   = NULL;
  double *FADIIS         = NULL;
  double *FBDIIS         = NULL;
  double *work           = NULL;

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
  RealMap B(BMem,lenCoeff,lenCoeff);
  

  //lapack variables for DIIS
  ErrorAlphaMem = BMem + lenB;
  if(!this->RHF_) {
    ErrorBetaMem = ErrorAlphaMem + lenF*(lenCoeff-1);
    FADIIS = ErrorBetaMem + lenF*(lenCoeff-1);
    FBDIIS = FADIIS + lenF*(lenCoeff-1);
    coef = FBDIIS + lenF*(lenCoeff-1);
  } else {
    FADIIS = ErrorAlphaMem + lenF*(lenCoeff-1);
    coef = FADIIS + lenF*(lenCoeff-1);
  }
  int *iPiv = new int[lenCoeff];
  int row=lenCoeff;
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
    RealMap ErrA(ErrorAlphaMem + (iter%(lenCoeff-1))*lenF,n,n);
    ErrA = (*this->fockA_) * (*this->densityA_) * (*this->aointegrals_->overlap_);
    ErrA -= (*this->aointegrals_->overlap_) * (*this->densityA_) * (*this->fockA_);
    memcpy(FADIIS+(iter%(lenCoeff-1))*lenF,this->fockA_->data(),lenF*sizeof(double));
    if(!this->RHF_){
      RealMap ErrB(ErrorBetaMem + (iter%(lenCoeff-1))*lenF,n,n);
      ErrB = (*this->fockB_) * (*this->densityB_) * (*this->aointegrals_->overlap_);
      ErrB -= (*this->aointegrals_->overlap_) * (*this->densityB_) * (*this->fockB_);
      memcpy(FBDIIS+(iter%(lenCoeff-1))*lenF,this->fockB_->data(),lenF*sizeof(double));
    }
    
    
    if(iter % (lenCoeff-1) == 5 && iter != 0){
      for(auto j = 0; j < (lenCoeff-1); j++)
      for(auto k = 0; k <= j          ; k++){
        RealMap EJA(ErrorAlphaMem + (j%(lenCoeff-1))*lenF,n,n);
        if(k==0) prettyPrint(this->fileio_->out,EJA,"Error "+std::to_string(j));
        RealMap EKA(ErrorAlphaMem + (k%(lenCoeff-1))*lenF,n,n);
        B(j,k) = -EJA.frobInner(EKA);
        if(!this->RHF_){
          RealMap EJB(ErrorBetaMem + (j%(lenCoeff-1))*lenF,n,n);
          RealMap EKB(ErrorBetaMem + (k%(lenCoeff-1))*lenF,n,n);
          B(j,k) = -EJB.frobInner(EKB);
        }
        B(k,j) = B(j,k);
      }
      for (auto l=0;l<lenCoeff-1;l++){
         B(lenCoeff-1,l)=-1.0;
	 B(l,lenCoeff-1)=-1.0;
      }
      B(lenCoeff-1,lenCoeff-1)=0;
      prettyPrint(this->fileio_->out,B,"B");
//    memset(coef,0.0,lenCoeff-1*sizeof(double)); // Zero out DIIS coef
      for(auto k = 0; k < lenCoeff;k++) coef[k] = 0.0; 
      coef[lenCoeff-1]=-1.0;
      this->fileio_->out << "COEFF Before" << endl;
      for(auto k = 0; k < lenCoeff; k++) this->fileio_->out << coef[k] << endl;
      dgesv_(&row,&nrhs,B.data(),&row, iPiv, coef,&row, &info);
      this->fileio_->out << "COEFF" << endl;
      for(auto k = 0; k < lenCoeff; k++) this->fileio_->out << coef[k] << endl;
      this->fockA_->setZero();
      if(!this->RHF_) this->fockB_->setZero();
      for(auto j = 0; j < lenCoeff-1; j++) {
        RealMap FA(FADIIS + (j%(lenCoeff-1))*lenF,n,n);
        *this->fockA_ += coef[j]*FA;
        prettyPrint(this->fileio_->out,FA,"Fock "+std::to_string(j));
        if(!this->RHF_) {
          RealMap FB(FBDIIS + (j%(lenCoeff-1))*lenF,n,n);
          *this->fockB_ += coef[j]*FB;
        }
      }
      prettyPrint(this->fileio_->out,*this->fockA_,"Total Fock");
    }

    PAlphaRMS=((*this->densityA_)-POldAlpha).norm();
    if(!this->RHF_) PBetaRMS = ((*this->densityB_) - POldBeta).norm();
    EDelta= this->totalEnergy-EOld;
    if(this->RHF_) this->printDensityInfo(PAlphaRMS,EDelta);     
    else this->printDensityInfo(PAlphaRMS,PBetaRMS,EDelta);     
    
    if(this->RHF_){
      if(PAlphaRMS < Dtol && pow((this->totalEnergy-EOld),2)<Etol) break;
    } else {
      if(PAlphaRMS < Dtol && PBetaRMS < Dtol && pow((this->totalEnergy-EOld),2)<Etol) break;
    }
  };
  
  //freeing the memory
  delete [] SCR;
  delete [] iPiv;
  
  if(iter >= maxIter)
    this->fileio_->out << "SCF Failed to converge within maximum number of iterations" << endl << endl;
  this->fileio_->out <<"\n"<<endl; 
  this->fileio_->out << bannerEnd <<endl;
  this->fileio_->out << "\nRequested convergence on RMS density matrix = " <<std::setw(5)<<Dtol <<"  within  " <<maxIter <<"  cycles."<<endl;
  this->fileio_->out << "Requested convergence on             energy = " <<Etol << endl;
  if(maxIter > iter){
    if(this->RHF_) this->fileio_->out << "\nSCF Done: E(RHF) = "<< this->totalEnergy << "  Eh  after  "<< iter+1 << "  cycles" <<endl;
    else this->fileio_->out << "\nSCF Done: E(UHF) = "<< this->totalEnergy << "  Eh  after  "<< iter+1 << "  cycles" <<endl;
  }
  this->fileio_->out << bannerEnd <<endl;
}; 
} // namespace ChronusQ
