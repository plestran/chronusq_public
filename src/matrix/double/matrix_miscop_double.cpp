/* 
 *  The Chronus Quantum (ChronusQ) software package is high-performace 
 *  computational chemistry software with a strong emphasis on explictly 
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
#include "matrix.h"
using ChronusQ::Matrix;
namespace ChronusQ {
template<> void Matrix<double>::pack(){
  double *tmp = new double[this->rows_*(this->rows_+1)/2];
  int iop = 1;
  char uplo = 'L';
  pkconv_(&iop,&uplo,&this->rows_,this->data_,&this->rows_,tmp);
  this->format_ = 1;
  this->len_ = this->rows_*(this->rows_+1)/2;
  delete this->data_;
  this->data_ = tmp;
};
template<> void Matrix<double>::unpack(){
  double *tmp = new double[this->rows_*this->cols_];
  int iop = 2;
  char uplo = 'L';
  pkconv_(&iop,&uplo,&this->rows_,tmp,&this->rows_,this->data_);
  this->format_ = 0;
  this->len_ = this->rows_*this->rows_;
  delete this->data_;
  this->data_ = tmp;
};
//sajan
template<> void Matrix<double>::norm_col(){
  double intermediate, currentData;
  double sum ;
  for (int i=0; i<this->cols_;i++) {
    sum = 0;
    for (int j=0; j<this->rows_;j++) {
      intermediate = data_[i*this->rows_+j];
      sum += intermediate*intermediate;
    };
    sum = (double)sqrt(sum);
    for (int k=0; k<this->rows_;k++) {
      currentData=data_[i*this->rows_+k];
      data_[i*this->rows_+k]=currentData/sum;
    };
  };
};
template<>
double* Matrix<double>::eigenvalue(){                // Return the default eigenvalue array
  if(haveEigen_==0) throw 3000;
  if(this->symm_=='G') return this->eigenvalue_re_;
  if(this->symm_=='S') return this->eigenvalue_;
};
template<>
double* Matrix<double>::eigenvalue_re(){
  if(haveEigen_==0) throw 3000;
  if(this->symm_=='G') return this->eigenvalue_re_;
  if(this->symm_=='S') return this->eigenvalue_;
}
template<>
double* Matrix<double>::eigenvalue_im(){
  double * zeros;
  if(haveEigen_==0) throw 3000;
  if(!this->realCmplxEig_){
    zeros = new double[this->rows_]; memset(zeros,0,len_*sizeof(double));
    return zeros;
  }else if(this->symm_=='G') return this->eigenvalue_im_;
  
}
/*
double* Matrix<double>::eigenvector_l(){
  if(this->haveEigen_!=1){ throw 3009;}
  else if(this->symm_=='G'){ return this->eigenvector_l_;}
  else if(this->symm_=='S'){ return this->eigenvector_;};
};
double* Matrix<double>::eigenvector_r(){
  if(this->haveEigen_!=1){ throw 3010;}
  else if(this->symm_=='G'){ return this->eigenvector_r_;}
  else if(this->symm_=='S'){ return this->eigenvector_;};
};
double* Matrix<double>::eigenvector(){
  if(this->haveEigen_!=1){ throw 3018;}
  else if(this->symm_=='S'){ return this->eigenvector_;}
  else if(this->symm_=='G'){ return this->eigenvector_r_;};
};
double* Matrix<double>::eigenvalue(){
  if(this->haveEigen_!=1){ throw 3019;}
  else if(this->symm_=='S'){ return this->eigenvalue_;}
  else if(this->symm_=='G'){ return this->eigenvalue_re_;};
};
*/
template<> void Matrix<double>::scaleDag(double x) {
  for(int i=0;i<cols_;i++) {
    data_[i*this->rows_ +i] = x*data_[i*this->rows_ +i];
  }
};

//
// Uses old fortran code to sort eigenvalues
//
template<> void Matrix<double>::eSort(){
  int iop = 1;

  if(this->symm_=='G') { 
    eigsrt_(&iop,this->eigenvalue_re_,this->eigenvalue_im_,this->eigenvector_r_,this->eigenvector_l_,&this->rows_,&this->rows_);
  } else if(this->symm_=='S') {
    eigsrt_(&iop,this->eigenvalue_,this->eigenvalue_,this->eigenvector_,this->eigenvector_,&this->rows_,&this->rows_);
  };
};

/* Attempt to generalize EigSrt to C++ (not tested, need to test)
template<> void Matrix<double>::eSort() {
  double tmp;
  double *W, *VR, *VL,
  if(this->symm_=='G') {
    W = this->eigenvalue_re_;
    VR = this->eigenvector_r_;
    VL = this->eigenvector_l_;
  } else if(this->symm_=='S') {
    W = this->eigenvalue_;
    VR = this->eigenvector_;
  }
  for(int i = 0; i < this->rows_-1; i++) {
    for(int j = i+1; j < this->rows_; j++) {
      if(W[i] > W[j]) {
        tmp = W[i];
        W[i] = W[j];
        W[j] = tmp;
        for(int k = 0; k < this->rows_; k++) {
          tmp = VR[k + i*this->rows_];
          VR[k + i*this->rows_] = VR[k + j*this->rows_];
          VR[k + j*this->rows_] = tmp;
        } // for k
      } // endIf
    } // for j
  } // for i
}
*/

template<> void Matrix<double>::allocEigen(){
  this->cleanEigen();
  if(this->symm_=='G'){

    // Allocate space for Eigensystem
    if(this->JOBVL_=='V') this->eigenvector_l_ = new (nothrow) double[this->len_];
    if(this->JOBVR_=='V') this->eigenvector_r_ = new (nothrow) double[this->len_];
    this->eigenvalue_re_ = new (nothrow) double[this->cols_];
    this->eigenvalue_im_ = new (nothrow) double[this->cols_];

    // Check if space was allocated without errors    
    if(this->eigenvalue_re_==NULL || this->eigenvalue_im_==NULL) { throw 3000;};
    if(this->eigenvector_l_==NULL || this->eigenvector_r_==NULL) { throw 3000;};

  } else if(this->symm_=='S') {

    // Allocate space for Eigensystem
    if(this->JOBVL_=='V' || this->JOBVR_=='V') this->eigenvector_ = new (nothrow) double[this->len_];
    this->eigenvalue_ = new (nothrow) double[this->cols_];

    // Check if space was allocated without errors    
    if(this->eigenvalue_ ==NULL) { throw 3000;};
    if(this->eigenvector_==NULL) { throw 3000;};
  } else throw 3000;
};



template<> int Matrix<double>::getLWORK(bool doGEP){
  double * test = new double[1];
  int LWORK = -1;
  int INFO = 0;
  char UPLO = 'U';
  
  if(doGEP){ throw 3000;
  } else {
    if(this->symm_=='G'){
      dgeev_(&this->JOBVL_,&this->JOBVR_,&this->rows_,this->data_,&this->rows_,this->eigenvalue_re_,
           this->eigenvalue_im_,this->eigenvector_l_,&this->rows_,this->eigenvector_r_,
           &this->rows_,test,&LWORK,&INFO);
    } else if(this->symm_=='S') {
      dsyev_(&this->JOBVR_,&UPLO,&this->rows_,this->eigenvector_,&this->rows_,this->eigenvalue_,
         test,&LWORK,&INFO);
    } else throw 3000;
  };
  LWORK = (int)test[0];
  return LWORK;
}
template<>  double** Matrix<double>::allocEigScr(bool doGEP, int &LWORK, Matrix *B){
  LWORK = this->getLWORK(doGEP);
  double** SCR = new double*[3];
  SCR[0] = NULL;
  SCR[1] = NULL;
  SCR[2] = NULL;

  SCR[0] = new (nothrow) double[LWORK];
  if(SCR[0]==NULL) throw 3000;
  if(this->symm_=='G') {
    SCR[1] = new double[this->len_];
    if(SCR[1]==NULL) throw 3000;
    for(int i=0;i<this->len_;i++) SCR[1][i] = this->data_[i];
  } else SCR[1] == NULL;

  if(this->symm_=='S'){
    for(int i=0;i<this->len_;i++) this->eigenvector_[i] = this->data_[i];
  };

  if(doGEP){
    SCR[2] = new double[this->len_];
    if(SCR[2]==NULL) throw 3000;
    for(int i=0;i<this->len_;i++) SCR[2][i] = B->data_[i];
  } 
  return SCR;
};
} // namespace ChronusQ
