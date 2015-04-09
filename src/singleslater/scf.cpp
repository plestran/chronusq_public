#include "singleslater.h"
using ChronusQ::Matrix;
//----------------------------------------//
// do the SCF                             //
// Sajan                                  //
//--------------- ------------------------//
void SingleSlater::SCF(){
  int n=this->nBasis_; //basis of the overlap matrix
  
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
  this->aointegrals_->overlap_->printAll();

  (*this->aointegrals_->overlap_)^2; 
  

  //(*transformX_) = (*diagOverlap_)^-0.5; 
  //transformX_->printAll();

//this->aointegrals_->overlap_->printAll();

//(*transformfockA_) = this->fockA_->transTN(*transformX_);

  
  
  /*
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

