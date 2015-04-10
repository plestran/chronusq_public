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
using ChronusQ::Matrix;
void Matrix<double>::cleanEigen(){
  delete[] this->eigenvalue_;
  delete[] this->eigenvector_;
  delete[] this->eigenvalue_re_;
  delete[] this->eigenvalue_im_;
  delete[] this->eigenvector_r_;
  delete[] this->eigenvector_l_;
  delete[] this->eigenvaluez_;
}
void Matrix<double>::allocEigen(){
  this->eigenvalue_ = new double[this->rows_];
  if(this->JOBVR_=='V') this->eigenvector_ = new double[this->len_];
  if(this->eigenvalue_==NULL)  throw 3100;
  if(this->eigenvector_==NULL && this->JOBVR=='V') throw 3101;
  if(this->symm_=='G'){
    this->eigenvalue_im_ = new double[this->rows_];
    if(this->JOBVL=='V') this->eigenvecltor_l_ = new double[this->len_];
    if(this->eigenvalue_im_==NULL)  throw 3102;
    if(this->eigenvector_l_==NULL && this->JOBVL=='V')  throw 3103;
  };

  this->eigenvector_r_ = this->eigenvector_;
  this->eigenvalue_re_ = this->eigenvalue_;
  if(this->symm_=='S')   this->eigenvector_l_ = this->eigenvector_;
}
double *Matrix<double>::allocScratch(){
  double* CPY = new double[this->len_];
  for(int i=0;i<this->len_;i++) CPY[i] = this->data_[i];
  return CPY;
}
int Matrix<double>::eigLWORK(int code, char UPLO){
  double *WORK = new double[1];
  int INFO,LWORK;
  if(     this->symm_=='G' && code==1){ 
    dgeev_(&this->JOBVL_,this->JOBVR_,&this->rows_,this->data_,&this->rows_,
      this->eigenvalue_re_,this->eigenvalue_im,this->eigenvector_l_;&this->rows_,
      this->eigenvector_r_,&this->rows_,WORK,-1,&INFO);
  }
  else if(this->symm_=='S' && code==1){ 
    dsyev_(&this->JOBVR_,&UPLO,&this->rows_,this->data_,&this->rows_,
      this->eigenvalue_,WORK,-1,&INFO);
  };
  LWORK = (int)WORK[0];
  return LWORK;
}
