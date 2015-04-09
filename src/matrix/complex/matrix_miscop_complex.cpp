#include "matrix.h"
using ChronusQ::Matrix;

namespace ChronusQ {
/*
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
*/
template<>
dcomplex* Matrix<dcomplex>::eigenvalue(){                // Return the default eigenvalue array
  if(haveEigen_==0) throw 3000;
  return this->eigenvalue_;
};
template<>
double* Matrix<dcomplex>::eigenvalue_re(){
  if(haveEigen_==0) throw 3000;
  return this->eigenvalue_re_;
}
template<>
double* Matrix<dcomplex>::eigenvalue_im(){
  if(haveEigen_==0) throw 3000;
  return this->eigenvalue_im_;
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
/*
template<> void Matrix<double>::scaleDag(double x) {
  for(int i=0;i<cols_;i++) {
    data_[i*this->rows_ +i] = x*data_[i*this->rows_ +i];
  }
};
template<> void Matrix<double>::eSort(){
  int iop = 1;

  if(this->symm_=='G') { 
    eigsrt_(&iop,this->eigenvalue_re_,this->eigenvalue_im_,this->eigenvector_r_,this->eigenvector_l_,&this->rows_,&this->rows_);
  } else if(this->symm_=='S') {
    eigsrt_(&iop,this->eigenvalue_,this->eigenvalue_,this->eigenvector_,this->eigenvector_,&this->rows_,&this->rows_);
  };
};
*/
template<> void Matrix<dcomplex>::allocEigen(){
  this->cleanEigen();
  if(this->symm_=='G'){

    // Allocate space for Eigensystem
    if(this->JOBVL_=='V') this->eigenvector_l_ = new (nothrow) dcomplex[this->len_];
    if(this->JOBVR_=='V') this->eigenvector_r_ = new (nothrow) dcomplex[this->len_];
    this->eigenvalue_ = new (nothrow) dcomplex[this->cols_];

    // Check if space was allocated without errors    
    if(this->eigenvector_l_==NULL || this->eigenvector_r_==NULL) { throw 3000;};
    if(this->eigenvalue_==NULL) { throw 3000;};

  } else if(this->symm_=='H') {

    // Allocate space for Eigensystem
    if(this->JOBVL_=='V' || this->JOBVR_=='V') this->eigenvector_ = new (nothrow) dcomplex[this->len_];
    this->eigenvalue_re_ = new (nothrow) double[this->cols_];

    // Check if space was allocated without errors    
    if(this->eigenvalue_re_ ==NULL) { throw 3000;};
    if(this->eigenvector_==NULL) { throw 3000;};
  } else throw 3000;
};



template<> int Matrix<dcomplex>::getLWORK(bool doGEP){
  dcomplex * test = new dcomplex[1];
  double *RWORK = new double[2*this->rows_];
  int LWORK = -1;
  int INFO = 0;
  char UPLO = 'U';
  
  if(doGEP){ throw 3000;
  } else {
    if(this->symm_=='G'){
      zgeev_(&this->JOBVL_,&this->JOBVR_,&this->rows_,this->data_,&this->rows_,this->eigenvalue_,
           this->eigenvector_l_,&this->rows_,this->eigenvector_r_, &this->rows_,test,&LWORK,RWORK,
           &INFO);
    } else if(this->symm_=='H') {
      zheev_(&this->JOBVR_,&UPLO,&this->rows_,this->eigenvector_,&this->rows_,this->eigenvalue_re_,
         test,&LWORK,RWORK,&INFO);
    } else throw 3000;
  };
  LWORK = (int)real(test[0]);
  delete[] test; delete[] RWORK;
  return LWORK;
}
template<>  dcomplex** Matrix<dcomplex>::allocEigScr(bool doGEP, int &LWORK, Matrix *B){
  LWORK = this->getLWORK(doGEP);
  dcomplex** SCR = new dcomplex*[3];
  SCR[0] = NULL;
  SCR[1] = NULL;
  SCR[2] = NULL;

  SCR[0] = new (nothrow) dcomplex[LWORK];
  if(SCR[0]==NULL) throw 3000;
  if(this->symm_=='G') {
    SCR[1] = new dcomplex[this->len_];
    if(SCR[1]==NULL) throw 3000;
    for(int i=0;i<this->len_;i++) SCR[1][i] = this->data_[i];
  } else SCR[1] == NULL;

  if(this->symm_=='H'){
    for(int i=0;i<this->len_;i++) this->eigenvector_[i] = this->data_[i];
  };

  if(doGEP){
    SCR[2] = new dcomplex[this->len_];
    if(SCR[2]==NULL) throw 3000;
    for(int i=0;i<this->len_;i++) SCR[2][i] = B->data_[i];
  } 
  cout << SCR << endl; 
  cout << SCR[0] << endl; 
  cout << SCR[1] << endl; 
  cout << SCR[2] << endl; 
  cout << " END " << endl;
  return SCR;
};
} // namespace ChronusQ
