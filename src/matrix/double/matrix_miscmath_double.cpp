#include "matrix.h"
/**************
 *  Addition  *
 **************/
template<> void Matrix<double>::add(const Matrix *a, const Matrix *b) {
  if(a==NULL||b==NULL) throw 3000;
  else if((this->len_!=a->len_)||(this->len_!=b->len_)||(a->len_!=b->len_)) throw 3000;

  for(int i=0;i<this->len_;i++) this->data_[i]=a->data_[i]+b->data_[i];
};
/*****************
 *  Subtraction  *
 *****************/
template<> void Matrix<double>::sub(const Matrix *a, const Matrix *b) {
  if(a==NULL||b==NULL) throw 3000;
  else if((this->len_!=a->len_)||(this->len_!=b->len_)||(a->len_!=b->len_)) throw 3000;

  for(int i=0;i<this->len_;i++) this->data_[i]=a->data_[i]-b->data_[i];
};
/*********************
 *  Diagonalization  *
 *********************/
template<> void Matrix<double>::doDiag(Matrix *B){
  bool doGEP = (B!=NULL);
  char UPLO = 'L';
  int LWORK,INFO;
  double **SCR;

  //  *Idiot Checks*  //
  // If Matrix is already diagonalized, early return
  if(this->haveEigen_!=0) return;
  // Diagonalization only works for Square matricies
  if(this->cols_!=this->rows_) throw 3000;
  // No packed Storage
  if(this->format_==1) throw 3000;

  // Allocate space for eigensystem and scratch
  this->allocEigen();
  SCR = this->allocEigScr(doGEP,LWORK,B);


  if(doGEP) throw 3000;
  else {
    if(this->symm_=='G'){
      dgeev_(&this->JOBVL_,&this->JOBVR_,&this->rows_,SCR[1],&this->rows_,this->eigenvalue_re_,
           this->eigenvalue_im_,this->eigenvector_l_,&this->rows_,this->eigenvector_r_,
           &this->rows_,SCR[0],&LWORK,&INFO);
    } else if(this->symm_=='S'){
      dsyev_(&this->JOBVR_,&UPLO,&this->rows_,this->eigenvector_,&this->rows_,this->eigenvalue_,
         SCR[0],&LWORK,&INFO);
      if(this->JOBVR_!='V'){ delete this->eigenvector_;};
    };
  };
  
  for(int i = 0; i < 3; i++) delete[] SCR[i];
  delete[] SCR;
  this->haveEigen_=1;
};
template<> void Matrix<double>::diag(Matrix *B){ this->doDiag(B);};
/*
void Matrix<double>::diag() {
  char UPLO = 'L';
  int LWORK;
  double *CPY;
  // If the matrix has already been diagonalized, state so and move on
  // (no error throw)
  if(this->haveEigen_==1){
     cout << "Matrix Object "<< this << " has already been diagonalized" << endl;
     return;
  };
  // Diagonalization only works for Square matricies
  if(this->cols_!=this->rows_) { throw 3011;}
  if(this->symm_=='G') { 
    if(this->eigenvector_l_==NULL && this->JOBVL_ == 'V') this->eigenvector_l_ = new (nothrow) double[this->len_];
    if(this->eigenvector_r_==NULL && this->JOBVR_ == 'V') this->eigenvector_r_ = new (nothrow) double[this->len_];
    // Having already made the decision to diagonalize, delete the eigenvalues if they already exist.
    if(this->eigenvalue_re_!=NULL) { delete[] this->eigenvalue_re_;};
    if(this->eigenvalue_im_!=NULL) { delete[] this->eigenvalue_im_;};
    // Allocate space for eigenvalues
    this->eigenvalue_re_ = new (nothrow) double[this->cols_];
    this->eigenvalue_im_ = new (nothrow) double[this->cols_];
    // If there were problems in the allocation for either the eigenvalues or eigenvectors, throw an error
    if(this->eigenvalue_re_==NULL || this->eigenvalue_im_==NULL) { throw 3013;};
    if(this->eigenvector_l_==NULL || this->eigenvector_r_==NULL) { throw 3014;};
    CPY  = new (nothrow) double[this->len_];     // DGEEV overwrites the matrix, so allocate space for a copy
    if(CPY== NULL){ throw 3016;};                // Throw an error if you can't allocate the space
    for(int i=0;i<len_;i++) CPY[i] = this->data_[i];
    // Define workspace for general matrix diagonalization (this should eventaully be done by a routine call
    // to determine optimum workspace as opposed to a "one-size-fits-all" LWORK = N**3
    LWORK = 4*this->cols_;
  } else if(this->symm_=='S' && this->format_==0) {
    // Check if the eigenvectors exist, if not, allocate them and set flags to calculate
    // the eigenvectors
    if(this->eigenvector_==NULL && this->JOBVR_ == 'V') this->eigenvector_ = new (nothrow) double[this->len_];
    // Having already made the decision to diagonalize, delete the eigenvalues if they already exist.
    if(this->eigenvalue_!=NULL) { delete[] this->eigenvalue_;};
    // Allocate space for eigenvalues
    this->eigenvalue_ = new (nothrow) double[this->cols_];
    // If there were problems in the allocation for either the eigenvalues or eigenvectors, throw an error
    if(this->eigenvalue_ ==NULL) { throw 3013;};
    if(this->eigenvector_==NULL) { throw 3014;};
    for(int i=0; i<this->len_; i++){           // DSYEV overwrites the matrix with the eigenvectors, so
      this->eigenvector_[i] = this->data_[i];  //  copy the matrix into the eigenvector array
    }
    // Define workspace for general matrix diagonalization (this should eventaully be done by a routine call
    // to determine optimum workspace as opposed to a "one-size-fits-all" LWORK = N**3
    LWORK = 3*this->cols_;
  } else {
    // Currently, only general/symmetric real matricies matricies are implemented for diagonalization.
    throw 3012;
  };
  
  if(this->format_!=0){ throw 3017;};          // Only works for full stored matricies
  
  double *WORK = new (nothrow) double[LWORK];  // Try to allocate work space
  if(WORK==NULL){ throw 3015;};                // If you can't, throw an error
  int INFO;                                    // Parameter to decide if LAPACK did its job.

  clock_t start,finish;
  start = clock();
  // This generalization is for eventual implementation of DSPEV
  if(this->format_==0){
    if(this->symm_=='G'){
//    cout << "Diagonalizing Matrix Object:"<<setw(5)<<" "<<" \""<<this->name_<<"\" at " << this << " using DGEEV" << endl;
      dgeev_(&this->JOBVL_,&this->JOBVR_,&this->rows_,CPY,&this->rows_,this->eigenvalue_re_,
           this->eigenvalue_im_,this->eigenvector_l_,&this->rows_,this->eigenvector_r_,
           &this->rows_,WORK,&LWORK,&INFO);
    } else if(this->symm_=='S'){
      cout << "Diagonalizing Matrix Object:"<<setw(5)<<" "<<" \""<<this->name_<<"\" at " << this << " using DSYEV" << endl;
      dsyev_(&this->JOBVR_,&UPLO,&this->rows_,this->eigenvector_,&this->rows_,this->eigenvalue_,
         WORK,&LWORK,&INFO);
      if(this->JOBVR_!='V'){ delete this->eigenvector_;};
    };
  };
  finish = clock();
//cout<<"Diagonalization of Matrix Object: \""<<this->name_<<"\" at " << this << " finished after "<<setprecision(8)<<(double)(finish-start)/CLOCKS_PER_SEC<<" CPU seconds."<<endl;
  delete[] WORK;
  if(this->symm_=='G'){ delete[] CPY;};
  this->haveEigen_=1;
};
void Matrix<double>::diag(Matrix *B){
  // Diagonalize "this" subject to the constraint B (i.e. Ax = wBx)
  char JOBVL = 'N';
  char JOBVR = 'N';
  char UPLO = 'L';
  int LWORK;
  int ITYPE = 1;
  double *CPYA, *CPYB, *ALPHAR, *ALPHAI, *BETA;
  // If the matrix has already been diagonalized, state so and move on
  // (no error throw)
  if(this->haveEigen_==1){
     cout << "Matrix Object "<< this << " has already been diagonalized" << endl;
     return;
  };
  // Diagonalization only works for Square matricies (both A and X)
  if((this->cols_!=this->rows_)||(B->cols_!=B->rows_)) { throw 3011;}
  if(this->symm_=='G') {
    // Check if the eigenvectors exist, if not, allocate them and set flags to calculate
    // the eigenvectors
    if(this->eigenvector_l_==NULL) {
      JOBVL = 'V';
      this->eigenvector_l_ = new (nothrow) double[this->len_];
    };
    if(this->eigenvector_r_==NULL) {
      JOBVR = 'V';
      this->eigenvector_r_ = new (nothrow) double[this->len_];
    };
    // Having already made the decision to diagonalize, delete the eigenvalues if they already exist.
    if(this->eigenvalue_re_!=NULL) { delete[] this->eigenvalue_re_;};
    if(this->eigenvalue_im_!=NULL) { delete[] this->eigenvalue_im_;};
    // Allocate space for eigenvalues
    this->eigenvalue_re_ = new (nothrow) double[this->cols_];
    this->eigenvalue_im_ = new (nothrow) double[this->cols_];
    // If there were problems in the allocation for either the eigenvalues or eigenvectors, throw an error
    if(this->eigenvalue_re_==NULL || this->eigenvalue_im_==NULL) { throw 3013;};
    if(this->eigenvector_l_==NULL || this->eigenvector_r_==NULL) { throw 3014;};
    CPYA  = new (nothrow) double[this->len_];    // DGEGV overwrites the A matrix, so allocate space for a copy
    CPYB  = new (nothrow) double[this->len_];    // DGEGV overwrites the B matrix, so allocate space for a copy
    if((CPYA== NULL)||(CPYB== NULL)){ throw 3016;};  // Throw an error if you can't allocate the space for copies
    for(int i=0;i<this->len_;i++){                     // Copy the matrix into the CPYA/B pointer
      CPYA[i] = this->data_[i];
      CPYB[i] = B->data_[i];
    };
    // Tempory storage for the possible overflow of eigenvalues by DGEGV. These will be combined to form the
    // real and imaginary parts of the eigenvalues in post processing.
    ALPHAR = new (nothrow) double[this->cols_];
    ALPHAI = new (nothrow) double[this->cols_];
    BETA   = new (nothrow) double[this->cols_];
    if(ALPHAR==NULL||ALPHAI==NULL||BETA==NULL){ throw 3018;};
    LWORK = 8*this->cols_;
  } else if(this->symm_=='S' && B->symm_=='S' && this->format_==0) {
    // Check if the eigenvectors exist, if not, allocate them and set flags to calculate
    // the eigenvectors
    if(this->eigenvector_==NULL) {
      JOBVR = 'V';
      this->eigenvector_ = new (nothrow) double[this->len_];
    };
    // Having already made the decision to diagonalize, delete the eigenvalues if they already exist.
    if(this->eigenvalue_!=NULL) { delete[] this->eigenvalue_;};
    // Allocate space for eigenvalues
    this->eigenvalue_ = new (nothrow) double[this->cols_];
    // If there were problems in the allocation for either the eigenvalues or eigenvectors, throw an error
    if(this->eigenvalue_ ==NULL) { throw 3013;};
    if(this->eigenvector_==NULL) { throw 3014;};
    CPYB  = new (nothrow) double[this->len_];    // DGEGV overwrites the B matrix, so allocate space for a copy
    if(CPYB== NULL){ throw 3016;};               // Throw an error if you can't allocate the space for copies
    for(int i=0; i<this->len_; i++){           // DSYEV overwrites the matrix with the eigenvectors, so
      this->eigenvector_[i] = this->data_[i];  //  copy the matrix into the eigenvector array
      CPYB[i] = B->data_[i];
    }
    // Define workspace for general matrix diagonalization (this should eventaully be done by a routine call
    // to determine optimum workspace as opposed to a "one-size-fits-all" LWORK = N**3
    LWORK = 3*this->cols_-1;
  } else {
    // Currently, only general/symmetric real matricies matricies are implemented for diagonalization.
    throw 3012;
  };
  if(this->format_!=0){ throw 3017;};          // Only works for full stored matricies
  
  double *WORK = new (nothrow) double[LWORK];  // Try to allocate work space
  if(WORK==NULL){ throw 3015;};                // If you can't, throw an error
  int INFO;                                    // Parameter to decide if LAPACK did its job.

  clock_t start,finish;
  start = clock();
  if(this->symm_=='G'){
    cout << "Diagonalizing Matrix Object:"<<setw(5)<<" "<<" \""<<this->name_<<"\" subject to " << B->name_ << " using DGEGV" << endl;
    dgegv_(&JOBVL,&JOBVR,&this->rows_,CPYA,&this->rows_,CPYB,&this->rows_,ALPHAR,ALPHAI,BETA,
         this->eigenvector_l_,&this->rows_,this->eigenvector_r_, &this->rows_,WORK,&LWORK,
         &INFO);
    for(int i=0;i<this->cols_;i++){ 
      this->eigenvalue_re_[i] = ALPHAR[i]/BETA[i];
      this->eigenvalue_im_[i] = ALPHAI[i]/BETA[i];
    };
  } else if(this->symm_=='S'){
    cout << "Diagonalizing Matrix Object:"<<setw(5)<<" "<<" \""<<this->name_<<"\" subject to " << B->name_ << " using DSYGV" << endl;
    dsygv_(&ITYPE,&JOBVR,&UPLO,&this->rows_,this->eigenvector_,&this->rows_,CPYB,&this->rows_,this->eigenvalue_,
       WORK,&LWORK,&INFO);
    
    if(JOBVR!='V'){ delete this->eigenvector_;};
  };
  finish = clock();
  cout<<"Diagonalization of Matrix Object: \""<<this->name_<<"\" at " << this << " finished after "<<setprecision(8)<<(double)(finish-start)/CLOCKS_PER_SEC<<" CPU seconds."<<endl;
  delete[] WORK; delete[] CPYB;
  if(this->symm_=='G'){ delete[] CPYA; delete[] ALPHAR; delete[] ALPHAI; delete[] BETA;};
  this->haveEigen_=1;
}
*/
/***************
 *  Transpose  *
 ***************/
template<>
void Matrix<double>::transposeHard(){
  double *tmp = new (nothrow) double[this->len_];
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
void Matrix<double>::adjointHard(){ this->transposeHard();};
//-----------------//
//      trace      //
//-----------------//
template<>
double Matrix<double>::trace() {
  double tmpVal = 0.0;
  if (this->rows_!=this->cols_) throw 3007;
  if (this->format_==0) {
    for(int i=0;i<this->rows_;i++) tmpVal+=this->data_[i*(this->rows_)+i];
  } else if (this->format_==1) {
    for(int i=0;i<this->rows_;i++) tmpVal+=this->data_[i*(this->rows_)-i*(i-1)/2];
  };
  return tmpVal;
};
template<>
double Matrix<double>::scalarProd(Matrix *m) {
  if(this->len_!=m->len_) throw 3006;
  double tmpVal = math.zero;
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
Matrix<double>::TNT Matrix<double>::transNT(const Matrix &X){
  TNT  s;
  s.a = this;
  s.x = &X;
  return s;
}
template<>
Matrix<double>::TTN Matrix<double>::transTN(const Matrix &X){
  TTN  s;
  s.a = this;
  s.x = &X;
  return s;
}
