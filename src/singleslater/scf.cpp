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
using ChronusQ::SingleSlater;
//----------------------------------------//
// do the SCF                             //
// Sajan                                  //
//--------------- ------------------------//
double E_delta;
double P_Rms;

void SingleSlater::SCF(){
  if(!this->aointegrals_->haveAOOneE) this->aointegrals_->computeAOOneE();
  double E_old;
  int maxIte    = 128; 
  int n=this->nBasis_; 
  double Dtol = 1e-10;
  double Etol = 1e-8;

  RealMatrix          X(n,n);
  RealMatrix  	     Fp(n,n);
  RealMatrix  	     Cp(n,n);
  RealMatrix      P_old(n,n);
  
  
  int i;
  X=(*this->aointegrals_->overlap_).pow(-0.5);
  for (i=0; i<maxIte;i++){
    this->fileio_->out << endl << endl << bannerTop <<endl;  
    this->fileio_->out << "SCF iteration:"<< i+1 <<endl;  
    this->fileio_->out << bannerEnd <<endl;  
    
    P_old = (*this->densityA_);
    E_old = this->totalEnergy;
    Fp    = X.transpose()*(*this->fockA_)*X;
    Eigen::SelfAdjointEigenSolver<RealMatrix> sys(Fp);
    Cp    = sys.eigenvectors();
    
    (*this->epsA_)= sys.eigenvalues();
    (*this->moA_) = X*Cp;
    this->formDensity();
    this->formFock();
    this->computeEnergy();
     
    P_Rms=((*this->densityA_)-P_old).norm();
    E_delta= this->totalEnergy-E_old;
    this->printDensityinf();     

    if(((*this->densityA_)-P_old).norm()<Dtol && pow((this->totalEnergy-E_old),2)<Etol){break;};
  };
  
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
