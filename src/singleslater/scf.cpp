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
#include <eiginterface.h>
using ChronusQ::Controls;
using ChronusQ::FileIO;
using ChronusQ::SingleSlater;

extern "C" void dgesv_(int *n, int *nrhs, double *A, int *lda , int *ipiv, double *b, int *ldb, int *info);
extern "C" void dsyev_(char *Jobz, char *Uplo, int *n, double *A, int *lda, double *W, double *work, int *lwork, int *info);
//----------------------------------------//
// do the SCF                             //
// Sajan                                  //
//----------------------------------------//
double E_delta;
double P_Rms;

void SingleSlater::SCF(){
  if(!this->aointegrals_->haveAOOneE) this->aointegrals_->computeAOOneE();
  double E_old;
  int maxIte    = 128; 
  int n=this->nBasis_; 
  double Dtol = 1e-10;
  double Etol = 1e-8;
  /*
  RealMatrix          X(n,n);
  RealMatrix  	     Fp(n,n);
  RealMatrix  	     Cp(n,n);
  RealMatrix      P_old(n,n);
  RealMatrix          B(7,7); 
  */
  std::vector<RealMatrix> Error;
  std::vector<RealMatrix> fock;
  int i; 
  //RealMatrix memory allocation
  int lenX = n*n; // X
  int lenFp = n*n; // Fp
  int lenPold = n*n; // Pold
  int lenB = 49;  // B
  int lenCoeff = 7;   // Coeff
  int lenEigV = n;   // Eigenvalues
  int lwork = 4*n; // LAPACK Workspace
  int LenScr = lenX + lenFp + lenPold + lenB + lenCoeff + lenEigV + lwork;

  double *SCR, *X_m, *Fp_m, *P_old_m, *B_m, *coef, *eig_values, *work;

  SCR = new double [LenScr];
  X_m= SCR;
  Fp_m=X_m + lenX;
  P_old_m= Fp_m + lenFp;
  B_m=P_old_m + lenPold;

  RealMap X(X_m,n,n);
  RealMap Fp(Fp_m,n,n);
  RealMap P_old(P_old_m,n,n);
  RealMap B(B_m,7,7);
  

  //lapack variables for DIIS
  coef = B_m + lenB;
  int *ipiv = new int[7];
  int row=7;
  int nrhs=1;
  int info=-1;

  //lapack variables for F'C=CE
  char j='V';
  char u='U';
  eig_values = coef + lenCoeff;
  work = eig_values + lenEigV;
   
  
  X=(*this->aointegrals_->overlap_).pow(-0.5);
  for (i=0; i<maxIte;i++){
    this->fileio_->out << endl << endl << bannerTop <<endl;  
    this->fileio_->out << "SCF iteration:"<< i+1 <<endl;  
    this->fileio_->out << bannerEnd <<endl;  
    
    P_old = (*this->densityA_);
    E_old = this->totalEnergy;
    Fp    = X.transpose()*(*this->fockA_)*X;
    dsyev_(&j,&u,&n, Fp.data(), &n, eig_values, work, &lwork, &info);
    Fp.transposeInPlace();
    (*this->moA_) = X*Fp;
    prettyPrint(this->fileio_->out,(*this->moA_),"Alpha Fock");
    this->formDensity();
    this->formFock();
    this->computeEnergy();
    
    
    //Implimenting DIIS
    if(i % 6 ==0){
      Error.push_back((*this->fockA_)*(*this->densityA_)*(*this->aointegrals_->overlap_)-(*this->aointegrals_->overlap_)*(*this->densityA_)*(*this->fockA_));
      fock.push_back(*this->fockA_);
    }
    if(i % 6 ==1){
      Error.push_back((*this->fockA_)*(*this->densityA_)*(*this->aointegrals_->overlap_)-(*this->aointegrals_->overlap_)*(*this->densityA_)*(*this->fockA_));
      fock.push_back(*this->fockA_);
    }
    if(i % 6 ==2){
      Error.push_back((*this->fockA_)*(*this->densityA_)*(*this->aointegrals_->overlap_)-(*this->aointegrals_->overlap_)*(*this->densityA_)*(*this->fockA_));
      fock.push_back(*this->fockA_);
    }
    if(i % 6 ==3){
      Error.push_back((*this->fockA_)*(*this->densityA_)*(*this->aointegrals_->overlap_)-(*this->aointegrals_->overlap_)*(*this->densityA_)*(*this->fockA_));
      fock.push_back(*this->fockA_);
    }
    if(i % 6 ==4){
      Error.push_back((*this->fockA_)*(*this->densityA_)*(*this->aointegrals_->overlap_)-(*this->aointegrals_->overlap_)*(*this->densityA_)*(*this->fockA_));
      fock.push_back(*this->fockA_);
    }
    if(i % 6 ==5){
      Error.push_back((*this->fockA_)*(*this->densityA_)*(*this->aointegrals_->overlap_)-(*this->aointegrals_->overlap_)*(*this->densityA_)*(*this->fockA_));
      fock.push_back(*this->fockA_);
    }
    
    
    if(i % 6==0 && i!=0){
      
      for (auto j=0;j<Error.size();j++){
        for (auto k=0; k<=j;k++){
          B(j,k)=(Error[j]*(Error[k].transpose())).trace();
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
      dgesv_(&row,&nrhs,B.data(),&row, ipiv, coef,&row, &info);
      RealMatrix interme(n,n);
      for (auto j=0;j<6;j++){
        interme = interme+ (coef[j]*fock[j]);
      }
      *this->fockA_=interme;
      Error.clear();
      fock.clear();
    }

    P_Rms=((*this->densityA_)-P_old).norm();
    E_delta= this->totalEnergy-E_old;
    this->printDensityinf();     
     
    if(((*this->densityA_)-P_old).norm()<Dtol && pow((this->totalEnergy-E_old),2)<Etol){break;};
  };
  
  //freeing the memory
  /*
  delete[] coef;
  delete[] ipiv;
  delete[] X_m;
  delete[] P_old_m;
  delete[] Fp_m;
  delete[] work;
  delete[] B_m;
  delete[] eig_values;
  */
  this->fileio_->out <<"\n"<<endl; 
  this->fileio_->out << bannerEnd <<endl;
  this->fileio_->out << "\nRequested convergence on RMS density matrix = " <<std::setw(5)<<Dtol <<"  within  " <<maxIte <<"  cycles."<<endl;
  this->fileio_->out << "Requested convergence on             energy = " <<Etol << endl;
  this->fileio_->out << "\nSCF Done: E(RHF) = "<< this->totalEnergy << "  Eh  after  "<< i+1 << "  cycles" <<endl;
  this->fileio_->out << bannerEnd <<endl;
}; 

void SingleSlater::printDensityinf(){
  this->fileio_->out<<"\nSCF Information:"<<endl;
  this->fileio_->out<<std::right<<std::setw(30)<<"    Delta-E = "<<std::setw(15)<<E_delta<<std::setw(5)<<" Eh "<<endl;
  this->fileio_->out<<std::right<<std::setw(30)<<"RMS Density = "<<std::setw(15)<<P_Rms<<endl;
};

