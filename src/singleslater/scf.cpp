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
void SingleSlater::SCF(){
  int n=this->nBasis_; //basis of the overlap matrix
  
  /*
  //Declaring new pointers
  double oldEnergy;
  double diffEnergy=0.1;
  double rms_Density;
  Matrix<double> *oldDensityA_=new Matrix<double>(n,n);
  Matrix<double> *diffDensity=new Matrix<double>(n,n);
  Matrix<double> *diagOverlap_=new Matrix<double>(n,n);
  Matrix<double> *transformX_=new Matrix<double>(n,n);
  Matrix<double> *transformX_prime=new Matrix<double>(n,n);
  Matrix<double> *transformfockA_=new Matrix<double>(n,n);
  Matrix<double> *eig_fockA_=new Matrix<double>(n,n);
  Matrix<double> *rms_diffDensity=new Matrix<double>(n,n);
  
  //overlap matrix
  this->aointegrals_->overlap_->setSymm('S');
  this->aointegrals_->overlap_->unpack();

  (*this->aointegrals_->overlap_)^2; 
  

  //(*transformX_) = (*diagOverlap_)^-0.5; 
  //transformX_->printAll();

//this->aointegrals_->overlap_->printAll();

//(*transformfockA_) = this->fockA_->transTN(*transformX_);

  
  
  //calculate transformation matrix and its transpose
  (*transformX_)=(*diagOverlap_)^-0.5;
  transformX_->printAll();
  (*transformX_prime)=(*transformX_);
  transformX_prime->transpose(); 
  transformX_prime->printAll();
  
  for(int i=0;i<10;i++){
    fockA_->unpack();
    (*oldDensityA_) = (*this->densityA_);  	 	    //copy the old density to new one
    oldEnergy= this->totalEnergy;
    (*transformfockA_)=(*transformX_prime)*(*fockA_)
    (*transformfockA_)=(*transformfockA_)*(*transformX_);
    (*eig_fockA_)=transformfockA_->eigenvector();
    (*moA_)=(transformX_)*(*eig_fockA);		            //new MO coefficient matrix
    fockA_->pack();
    formDensity(fileio, controls);                          //form a new density matrix
    intTwoE(molecule,basis,this,fileio,controls);
    formFock(fileio, controls);                             //forming a new fock matrix
    computeEnergy(fileio, controls);                        //calcualte the new total energy with new density
    cout << this->totalEnergy << endl; 
    diffEnergy=((this->totalEnergy - oldEnergy)/oldEnergy)*100;
    (*diffDensity)=(*this->densityA_) - (*oldDensityA_);
    (*rms_diffDensity)=(*diffDensity)*(*diffDensity);
    //rms_Density=*(rms_diffDensity->sumAllelem());
    rms_Density=(double)sqrt(rms_Density);
  };

  */
};     

