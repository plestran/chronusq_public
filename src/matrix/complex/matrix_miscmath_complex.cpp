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
namespace ChronusQ {
/**************
 *  Addition  *
 **************/
template<> void Matrix<dcomplex>::add(const Matrix *a, const Matrix *b) {
  if(a==NULL||b==NULL) throw 3000;
  else if((this->len_!=a->len_)||(this->len_!=b->len_)||(a->len_!=b->len_)) throw 3000;

  for(int i=0;i<this->len_;i++) this->data_[i]=a->data_[i]+b->data_[i];
};
/*****************
 *  Subtraction  *
 *****************/
template<> void Matrix<dcomplex>::sub(const Matrix *a, const Matrix *b) {
  if(a==NULL||b==NULL) throw 3000;
  else if((this->len_!=a->len_)||(this->len_!=b->len_)||(a->len_!=b->len_)) throw 3000;

  for(int i=0;i<this->len_;i++) this->data_[i]=a->data_[i]-b->data_[i];
};
/*********************
 *  Diagonalization  *
 *********************/
template<> void Matrix<dcomplex>::doDiag(Matrix *B){
  bool doGEP = (B!=NULL);
  char UPLO = 'L';
  int LWORK,INFO;
  dcomplex **SCR;
  double *RWORK = NULL;

  //  *Idiot Checks*  //
  // If Matrix is already diagonalized, early return
  if(this->haveEigen_!=0) return;
  // Diagonalization only works for Square matricies
  if(this->cols_!=this->rows_) throw 3000;
  // No packed Storage
  if(this->format_==1) throw 3000;

  // Allocate space for eigensystem and scratch
  cout << " HERE 1" << endl;
  this->allocEigen();
  SCR = this->allocEigScr(doGEP,LWORK,B);
  RWORK = new double[2*this->rows_];
  cout << SCR << endl; 
  cout << SCR[0] << endl; 
  cout << SCR[1] << endl; 
  cout << SCR[2] << endl; 
  cout << RWORK << endl;
  cout << " HERE 2" << endl;


  if(doGEP) throw 3000;
  else {
    if(this->symm_=='G'){
      cout << " HERE 3" << endl;
      zgeev_(&this->JOBVL_,&this->JOBVR_,&this->rows_,SCR[1],&this->rows_,this->eigenvalue_,
           this->eigenvector_l_,&this->rows_,this->eigenvector_r_,&this->rows_,SCR[0],&LWORK,
           RWORK,&INFO);
    } else if(this->symm_=='H'){
      cout << SCR << endl; 
      cout << SCR[0] << endl; 
      cout << SCR[1] << endl; 
      cout << SCR[2] << endl; 
      cout << RWORK << endl;
      cout << " HERE 4" << endl;
      zheev_(&this->JOBVR_,&UPLO,&this->rows_,this->eigenvector_,&this->rows_,this->eigenvalue_re_,
         SCR[0],&LWORK,RWORK,&INFO);
      cout << SCR << endl; 
      cout << SCR[0] << endl; 
      cout << SCR[1] << endl; 
      cout << SCR[2] << endl; 
      cout << RWORK << endl;
      cout << " HERE 5" << endl;
      if(this->JOBVR_!='V'){ delete this->eigenvector_;};
    };
  };
  cout << SCR << endl; 
  cout << SCR[0] << endl; 
  cout << SCR[1] << endl; 
  cout << SCR[2] << endl; 
  cout << RWORK << endl;
  cout << " HERE65" << endl;
  for(int i = 0; i < 3; i++) delete[] SCR[i];
  cout << " HERE65" << endl;
  delete[] SCR; 
  delete[] RWORK;
  this->haveEigen_=1;
};
template<> void Matrix<dcomplex>::diag(Matrix *B){ this->doDiag(B);};
/***************
 *  Transpose  *
 ***************/
template<>
void Matrix<dcomplex>::transposeHard(){
  dcomplex *tmp = new (nothrow) dcomplex[this->len_];
  for(int i = 0; i < this->rows_; i++){
    for(int j = 0; j < this->cols_; j++){
      tmp[j*this->cols_ + i] = this->data_[i*this->rows_ + j];
    }
  }
  int swp = this->cols_;
  this->cols_ = this->rows_;
  this->rows_ = swp;
  delete[] data_;
  data_ = tmp;
};
template<>
void Matrix<dcomplex>::adjointHard(){
  dcomplex *tmp = new (nothrow) dcomplex[this->len_];
  for(int i = 0; i < this->rows_; i++){
    for(int j = 0; j < this->cols_; j++){
      tmp[j*this->cols_ + i] = conj(this->data_[i*this->rows_ + j]);
    }
  }
  int swp = this->cols_;
  this->cols_ = this->rows_;
  this->rows_ = swp;
  delete[] data_;
  data_ = tmp;
};
//-----------------//
//      trace      //
//-----------------//
template<>
dcomplex Matrix<dcomplex>::trace() {
  dcomplex tmpVal = 0.0;
  if (this->rows_!=this->cols_) throw 3007;
  if (this->format_==0) {
    for(int i=0;i<this->rows_;i++) tmpVal+=this->data_[i*(this->rows_)+i];
  } else if (this->format_==1) {
    for(int i=0;i<this->rows_;i++) tmpVal+=this->data_[i*(this->rows_)-i*(i-1)/2];
  };
  return tmpVal;
};
template<>
dcomplex Matrix<dcomplex>::scalarProd(Matrix *m) {
  if(this->len_!=m->len_) throw 3006;
  dcomplex tmpVal = math.zero;
  int i;
  if(this->format()==0) for(i=0;i<this->len();i++) tmpVal+=(this->data()[i])*(m->data()[i]);
  else if(this->format()==1){
    for(i=0;i<this->len();i++) {
      tmpVal+=math.two*(this->data()[i])*(m->data()[i]);
    };
    for(i=0;i<this->rows();i++) {
      tmpVal-=(*this)(i,i)*(*m)(i,i);
    };
  };
  return tmpVal;
};
template<>
Matrix<dcomplex>::TNT Matrix<dcomplex>::transNT(const Matrix &X){
  TNT  s;
  s.a = this;
  s.x = &X;
  return s;
}
template<>
Matrix<dcomplex>::TTN Matrix<dcomplex>::transTN(const Matrix &X){
  TTN  s;
  s.a = this;
  s.x = &X;
  return s;
}
} // namespace ChronusQ
