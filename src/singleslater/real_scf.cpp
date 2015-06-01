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
double E_delta;
double P_Rms;

namespace ChronusQ {
template<>
void SingleSlater<double>::printDensityinf(){
  this->fileio_->out<<"\nSCF Information:"<<endl;
  this->fileio_->out<<std::right<<std::setw(30)<<"    Delta-E = "<<std::setw(15)<<E_delta<<std::setw(5)<<" Eh "<<endl;
  this->fileio_->out<<std::right<<std::setw(30)<<"RMS Density = "<<std::setw(15)<<P_Rms<<endl;
};
template<>
void SingleSlater<double>::SCF(){
  if(!this->aointegrals_->haveAOOneE) this->aointegrals_->computeAOOneE();
  double EOld;
  int maxIter    = 250; 
  int n = this->nBasis_; 
  double Dtol = 1e-10;
  double Etol = 1e-8;
  std::vector<RealMatrix> ErrorAlpha;
  std::vector<RealMatrix> fockAlpha;
  int iter; 
  
  //RealMatrix memory allocation
  int lenX = n*n; // X
  int lenF = n*n; // FpAlpha
  int lenP = n*n; // Pold
  int lenB = 49;  // B
  int lenCoeff = 7;   // Coeff
  int lenEig = n;   // Eigenvalues
  int lWork = 4*n; // LAPACK Workspace
  int LenScr = lenX + lenF + lenP + lenB + lenCoeff + lenEig + lWork;
  if(!this->RHF_) LenScr += lenF + lenP + lenEig;

  double *SCR, *Xm, *FPmAlpha, *PoldmAlpha, *Bm, *coef, *eigValuesAlpha, *work;
  double *FPmBeta, *PoldmBeta, *eigValuesBeta;

  SCR = new double [LenScr]; // Allocated scratch space
  Xm         = SCR;
  FPmAlpha   = Xm + lenX;
  PoldmAlpha = FPmAlpha  +  (this->RHF_ + 1)*lenF;
  Bm         = PoldmAlpha + (this->RHF_ + 1)*lenP;
  if(this->RHF_){
    FPmBeta   = FPmAlpha   + lenF;
    PoldmBeta = PoldmAlpha + lenP;
  }

  RealMap X(Xm,n,n);
  RealMap FpAlpha(FPmAlpha,n,n);
  RealMap POldAlpha(PoldmAlpha,n,n);
  RealMap FpBeta(FPmBeta,n,n);
  RealMap POldBeta(PoldmBeta,n,n);
  RealMap B(Bm,7,7);
  

  //lapack variables for DIIS
  coef = Bm + lenB;
  int *iPiv = new int[7];
  int row=7;
  int nrhs=1;
  int info=-1;

  //lapack variables for F'C'=C'E
  char j='V';
  char u='U';
  eigValuesAlpha = coef + lenCoeff;
  work = eigValuesAlpha + (this->RHF_ + 1)*lenEig;
  if(this->RHF_) eigValuesBeta = eigValuesAlpha + lenEig;
   
  
  X=(*this->aointegrals_->overlap_).pow(-0.5);
  for (iter=0; iter<maxIter; iter++){
    this->fileio_->out << endl << endl << bannerTop <<endl;  
    this->fileio_->out << "SCF iteration:"<< iter+1 <<endl;  
    this->fileio_->out << bannerEnd <<endl;  
    
    POldAlpha = (*this->densityA_);
    if(!this->RHF_) POldBeta = (*this->densityA_);
    EOld = this->totalEnergy;
    FpAlpha    = X.transpose()*(*this->fockA_)*X;
    dsyev_(&j,&u,&n, FpAlpha.data(), &n, eigValuesAlpha, work, &lWork, &info);
    FpAlpha.transposeInPlace();
    (*this->moA_) = X*FpAlpha;
    prettyPrint(this->fileio_->out,(*this->moA_),"Alpha Fock");
    this->formDensity();
    this->formFock();
    this->computeEnergy();
    
    
    //Implimenting DIIS
    if(iter % 6 ==0){
      ErrorAlpha.push_back((*this->fockA_) * (*this->densityA_) * (*this->aointegrals_->overlap_) -
                      (*this->aointegrals_->overlap_) * (*this->densityA_) * (*this->fockA_));
      fockAlpha.push_back(*this->fockA_);
    }
    if(iter % 6 ==1){
      ErrorAlpha.push_back((*this->fockA_) * (*this->densityA_) * (*this->aointegrals_->overlap_) -
                      (*this->aointegrals_->overlap_) * (*this->densityA_) * (*this->fockA_));
      fockAlpha.push_back(*this->fockA_);
    }
    if(iter % 6 ==2){
      ErrorAlpha.push_back((*this->fockA_) * (*this->densityA_) * (*this->aointegrals_->overlap_) - 
                      (*this->aointegrals_->overlap_) * (*this->densityA_) * (*this->fockA_));
      fockAlpha.push_back(*this->fockA_);
    }
    if(iter % 6 ==3){
      ErrorAlpha.push_back((*this->fockA_) * (*this->densityA_) * (*this->aointegrals_->overlap_) - 
                      (*this->aointegrals_->overlap_) * (*this->densityA_) * (*this->fockA_));
      fockAlpha.push_back(*this->fockA_);
    }
    if(iter % 6 ==4){
      ErrorAlpha.push_back((*this->fockA_) * (*this->densityA_) * (*this->aointegrals_->overlap_) - 
                      (*this->aointegrals_->overlap_) * (*this->densityA_) * (*this->fockA_));
      fockAlpha.push_back(*this->fockA_);
    }
    if(iter % 6 ==5){
      ErrorAlpha.push_back((*this->fockA_) * (*this->densityA_) * (*this->aointegrals_->overlap_) - 
                      (*this->aointegrals_->overlap_) * (*this->densityA_) * (*this->fockA_));
      fockAlpha.push_back(*this->fockA_);
    }
    
    
    if(iter % 6==0 && iter!=0){
      
      for (auto j=0;j<ErrorAlpha.size();j++){
        cout << ErrorAlpha[j] << endl << endl;
        for (auto k=0; k<=j;k++){
          B(j,k)=(ErrorAlpha[j]*(ErrorAlpha[k].transpose())).trace();
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
      RealMatrix interme(n,n);
      for (auto j=0;j<6;j++){
        interme = interme+ (coef[j]*fockAlpha[j]);
      }
      *this->fockA_=interme;
      ErrorAlpha.clear();
      fockAlpha.clear();
    }

    P_Rms=((*this->densityA_)-POldAlpha).norm();
    E_delta= this->totalEnergy-EOld;
    this->printDensityinf();     
     
    if(((*this->densityA_)-POldAlpha).norm()<Dtol && pow((this->totalEnergy-EOld),2)<Etol){break;};
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
