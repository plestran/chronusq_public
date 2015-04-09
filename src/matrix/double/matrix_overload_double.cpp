#include "matrix.h"
using ChronusQ::Matrix;
namespace ChronusQ {
/**********
 *  Copy  *
 **********/
template<>
void Matrix<double>::operator=(const Matrix *m) {
  if(m==NULL){ throw 3000;}
  else if((this->cols_!=m->cols_)||(this->rows_!=m->rows_)){ throw 3000;}

  for(int i=0;i<this->len_;i++) this->data_[i]=m->data_[i];
};
template<>
void Matrix<double>::operator=(double *m) {
  if(m==NULL){ throw 3020;}

  for(int i = 0; i<this->len_;i++) this->data_[i]=m[i];
}
/**************
 *  Addition  *
 **************/
template<>
void Matrix<double>::operator+=(const Matrix *m) {
  if(m==NULL){ throw 3000; }
  else if((this->cols_!=m->cols_)||(this->rows_!=m->rows_)) { throw 3000; }

  for(int i=0;i<this->len_;i++) this->data_[i]+=m->data_[i];
};
/*****************
 *  Subtraction  *
 *****************/
template<>
void Matrix<double>::operator-=(const Matrix *m) {
  if(m==NULL){ throw 3000; }
  else if((this->cols_!=m->cols_)||(this->rows_!=m->rows_)){ throw 3000; }

  for(int i=0;i<this->len_;i++) this->data_[i]-=m->data_[i];
};
template<>
Matrix<double>::Diff Matrix<double>::operator-(const Matrix &m) const {
  Diff s;

  if(&m==NULL) throw 3000;
  else if((this->cols_!=(&m)->cols_)||(this->rows_!=(&m)->rows_)){ throw 3000; }

  s.a = this;
  s.b = &m;
  return s;
};
template<>
void Matrix<double>::operator=(const Diff& m){
  for(int i = 0; i < this->len_; i++){
    this->data_[i] = m.a->data_[i] - m.b->data_[i];
  };
};
/*************
 *  Product  *
 *************/
template<>
void Matrix<double>::operator=(const Product& m){
  if(m.a->cols_!=m.b->rows_){ throw 3008;};               // Incompatable dimensions
  if(m.a->format_!=0 && !m.a->vectorized_){ throw 3024;}; // Packed NYI

  double alpha = 1.0;
  double beta = 0.0;
  char UPLO = 'L';                                        // Assumes that 'S' is LT (TBR)

  //** Note that for symmetric matricies, only the lower triangle is referenced **//
  if(m.a->symm_=='G') {
    dgemm_(&m.a->trans_,&m.b->trans_,&m.a->rows_,&m.b->cols_,&m.a->cols_,&alpha,m.a->data_,
           &m.a->rows_,m.b->data_,&m.b->rows_,&beta,this->data_,&this->rows_);
  } else if(m.a->symm_=='S'){
    char SIDE = 'L';
    dsymm_(&SIDE,&UPLO,&m.a->rows_,&m.b->cols_,&alpha,m.a->data_,&m.a->rows_,m.b->data_,
           &m.b->rows_,&beta,this->data_,&this->rows_);
  } else if(m.b->symm_=='S'){
    char SIDE = 'R';
    dsymm_(&SIDE,&UPLO,&m.b->rows_,&m.a->cols_,&alpha,m.b->data_,&m.b->rows_,m.a->data_,
           &m.a->rows_,&beta,this->data_,&this->rows_);
  };
};
template<>
Matrix<double>::Product Matrix<double>::operator*(const Matrix &m) const{
  Product s;

  if(&m==NULL) throw 3000;

  s.a = this;
  s.b = &m;
  return s;
};
/**************************************
 *  Transformation Related Operators  *
 **************************************/ 
void doTNT(Matrix<double> *res, const Matrix<double> *a, const Matrix<double> *x) {
  if(a->cols_!=a->rows_){ throw 3025;};
  if(x->cols_!=a->rows_){ throw 3026;};
  if(a->format_!=0 || x->format_!=0){ throw 3027;};
  double *tmp = new double[a->len_];
  double alpha = 1.0; double beta = 0.0;
  char xtrans = 'N'; char tmptrans = 'N';
  if(a->symm_=='G') {
    dgemm_(&xtrans,&a->trans_,&x->rows_,&a->cols_,&x->cols_,&alpha,x->data_,
           &x->rows_,a->data_,&a->rows_,&beta,tmp,&x->rows_);
    xtrans = 'T';
    dgemm_(&tmptrans,&xtrans,&x->rows_,&a->cols_,&a->cols_,&alpha,tmp,
	 &x->rows_,x->data_,&x->rows_,&beta,res->data_,&res->rows_);
  } else if(a->symm_=='S'){
    char UPLO = 'L';
    char SIDE = 'R';
    dsymm_(&SIDE,&UPLO,&x->rows_,&a->cols_,&alpha,a->data_,&a->rows_,x->data_,
           &x->rows_,&beta,tmp,&x->rows_);
    xtrans = 'T';
    dgemm_(&tmptrans,&xtrans,&x->rows_,&x->rows_,&a->cols_,&alpha,tmp,
	 &x->rows_,x->data_,&x->rows_,&beta,res->data_,&res->rows_);
  };
  delete[] tmp;
};
template<> void Matrix<double>::operator=(const TNT& m){ doTNT(this,m.a,m.x);};

/* This is an old routine, but i'd like to keep it around just in case
template<> void Matrix<double>::operator=(const TNT& m){
  if(m.a->cols_!=m.a->rows_){ throw 3025;};
  if(m.x->cols_!=m.a->rows_){ throw 3026;};
  if(m.a->format_!=0 || m.x->format_!=0){ throw 3027;};
  double alpha = 1.0;
  double beta = 0.0;
  double *tmp = new double[m.a->len_];
  char xtrans = 'N';
  char tmptrans = 'N';
  if(m.a->symm_=='G') {
    dgemm_(&xtrans,&m.a->trans_,&m.x->rows_,&m.a->cols_,&m.x->cols_,&alpha,m.x->data_,
           &m.x->rows_,m.a->data_,&m.a->rows_,&beta,tmp,&m.x->rows_);
    xtrans = 'T';
    dgemm_(&tmptrans,&xtrans,&m.x->rows_,&m.a->cols_,&m.a->cols_,&alpha,tmp,
	 &m.x->rows_,m.x->data_,&m.x->rows_,&beta,this->data_,&this->rows_);
  } else if(m.a->symm_=='S'){
    char UPLO = 'L';
    char SIDE = 'R';
    dsymm_(&SIDE,&UPLO,&m.x->rows_,&m.a->cols_,&alpha,m.a->data_,&m.a->rows_,m.x->data_,
           &m.x->rows_,&beta,tmp,&m.x->rows_);
    xtrans = 'T';
    dgemm_(&tmptrans,&xtrans,&m.x->rows_,&m.x->rows_,&m.a->cols_,&alpha,tmp,
	 &m.x->rows_,m.x->data_,&m.x->rows_,&beta,this->data_,&this->rows_);
  };
  delete[] tmp;
};
*/
void doTTN(Matrix<double> *res, const Matrix<double> *a, const Matrix<double> *x) {
  if(a->cols_!=a->rows_){ throw 3028;};
  if(x->rows_!=a->rows_){ throw 3029;};
  if(a->format_!=0 || x->format_!=0){ throw 3030;};
  double *tmp = new double[a->len_];
  double alpha = 1.0; double beta = 0.0;
  char xtrans = 'T'; char tmptrans = 'N';
  if(a->symm_=='G') {
    dgemm_(&xtrans,&a->trans_,&x->cols_,&a->cols_,&x->rows_,&alpha,x->data_,
           &x->rows_,a->data_,&a->rows_,&beta,tmp,&x->cols_);
    xtrans = 'N';
    dgemm_(&tmptrans,&xtrans,&x->cols_,&a->cols_,&x->rows_,&alpha,tmp,
	 &x->cols_,x->data_,&x->rows_,&beta,res->data_,&res->rows_);
  } else if(a->symm_=='S'){
    char UPLO = 'L';
    char SIDE = 'L';
    // Because DSYMM doesn't take a transpose arguement for X, we have to compute
    // AX -> tmp. and then tell DGEMM to transpose (AX)T = (XT)(AT) = (XT)A
    // because AT=A
    dsymm_(&SIDE,&UPLO,&a->rows_,&x->cols_,&alpha,a->data_,&a->rows_,x->data_,
           &x->rows_,&beta,tmp,&a->rows_);
    tmptrans = 'T';
    xtrans = 'N';
    dgemm_(&tmptrans,&xtrans,&x->cols_,&x->cols_,&a->rows_,&alpha,tmp,
	 &a->rows_,x->data_,&x->rows_,&beta,res->data_,&res->rows_);
  };
  delete[] tmp;
};
template<> void Matrix<double>::operator=(const TTN& m){ doTTN(this,m.a,m.x);};
/* This is an old routine, but i'd like to keep it around just in case
template<> void Matrix<double>::operator=(const TTN& m){
  if(m.a->cols_!=m.a->rows_){ throw 3028;};
  if(m.x->rows_!=m.a->rows_){ throw 3029;};
  if(m.a->format_!=0 || m.x->format_!=0){ throw 3030;};
  double alpha = 1.0;
  double beta = 0.0;
  double *tmp = new double[m.a->len_];
  char xtrans = 'T';
  char tmptrans = 'N';
  if(m.a->symm_=='G') {
    dgemm_(&xtrans,&m.a->trans_,&m.x->cols_,&m.a->cols_,&m.x->rows_,&alpha,m.x->data_,
           &m.x->rows_,m.a->data_,&m.a->rows_,&beta,tmp,&m.x->cols_);
    xtrans = 'N';
    dgemm_(&tmptrans,&xtrans,&m.x->cols_,&m.a->cols_,&m.x->rows_,&alpha,tmp,
	 &m.x->cols_,m.x->data_,&m.x->rows_,&beta,this->data_,&this->rows_);
  } else if(m.a->symm_=='S'){
    char UPLO = 'L';
    char SIDE = 'L';
    // Because DSYMM doesn't take a transpose arguement for X, we have to compute
    // AX -> tmp. and then tell DGEMM to transpose (AX)T = (XT)(AT) = (XT)A
    // because AT=A
    dsymm_(&SIDE,&UPLO,&m.a->rows_,&m.x->cols_,&alpha,m.a->data_,&m.a->rows_,m.x->data_,
           &m.x->rows_,&beta,tmp,&m.a->rows_);
    tmptrans = 'T';
    xtrans = 'N';
    dgemm_(&tmptrans,&xtrans,&m.x->cols_,&m.x->cols_,&m.a->rows_,&alpha,tmp,
	 &m.a->rows_,m.x->data_,&m.x->rows_,&beta,this->data_,&this->rows_);
  };
  delete[] tmp;
};
*/
/*********
 *  A^x  *
 *********/
template<>
Matrix<double>::PWR Matrix<double>::operator^(const double &x) const{
  PWR s;
  s.a = this;
  s.x = &x;
  return s;
};
template<>
void Matrix<double>::operator=(const PWR &m){
  Matrix<double> *Vec =  new Matrix<double>(m.a->rows_,m.a->cols_);
  Matrix<double> *Val =  new Matrix<double>(m.a->rows_,m.a->cols_);
  double *tmp = new double[m.a->rows_];
  Vec->clearAll();
  Val->clearAll();
  if(m.a->symm_=='G')      for(int i=0; i<m.a->rows_; i++) tmp[i] = (double)pow(m.a->eigenvalue_re_[i],*m.x);
  else if(m.a->symm_=='S') for(int i=0; i<m.a->rows_; i++) tmp[i] = (double)pow(m.a->eigenvalue_[i],*m.x);
  Val->setDag(tmp);
  if(m.a->symm_=='G')      (*Vec) = m.a->eigenvector_r_;
  else if(m.a->symm_=='S') (*Vec) = m.a->eigenvector_;
  delete[] tmp;

  doTNT(this,Val,Vec);

  delete Vec; delete Val;
};
/* An old function, but I'd like to keep it around
void Matrix<double>::operator=(const PWR& m){
  Matrix *Vec =  new Matrix(m.a->rows_,m.a->cols_);
  Matrix *Val =  new Matrix(m.a->rows_,m.a->cols_);
  double *tmp =  new double[m.a->rows_];
  double *tmp2 = new double[m.a->len_]; 
  Vec->clearAll();
  Val->clearAll();
//if(m.a->haveEigen_==0) { m.a->diag()};
  for(int i=0; i<m.a->rows_; i++){ tmp[i] = (double)pow(m.a->eigenvalue_[i],*m.x);};
  (*Vec) = (m.a->eigenvector_);
  Val->setDag(tmp);
  Vec->printAll();
  Val->printAll();
  double alpha = 1.0; double beta = 0.0;
  char xtrans = 'N'; char tmptrans = 'N';
  dgemm_(&xtrans,&Val->trans_,&Vec->rows_,&Val->cols_,&Vec->cols_,&alpha,Vec->data_,
         &Vec->rows_,Val->data_,&Vec->cols_,&beta,tmp2,&Vec->rows_);
  xtrans = 'T';
  dgemm_(&tmptrans,&xtrans,&Val->rows_,&Val->cols_,&Val->cols_,&alpha,tmp2,
	 &Val->rows_,Vec->data_,&Val->cols_,&beta,this->data_,&Val->rows_);
  delete[] tmp; delete[] tmp2;
  delete Vec; delete Val;
};
*/
/************
 *  Exp(A)  *
 ************/
Matrix<double>::EXP exp(const Matrix<double>& m){
  Matrix<double>::EXP s;
  s.a = &m;
  return s;
};
template<>
void Matrix<double>::operator=(const EXP &m){
  Matrix<double> *Vec =  new Matrix<double>(m.a->rows_,m.a->cols_);
  Matrix<double> *Val =  new Matrix<double>(m.a->rows_,m.a->cols_);
  double *tmp = new double[m.a->rows_];
  Vec->clearAll();
  Val->clearAll();
  if(m.a->symm_=='G')      for(int i=0; i<m.a->rows_; i++) tmp[i] = (double)std::exp(m.a->eigenvalue_re_[i]);
  else if(m.a->symm_=='S') for(int i=0; i<m.a->rows_; i++) tmp[i] = (double)std::exp(m.a->eigenvalue_[i]);
  Val->setDag(tmp);
  Val->printAll();
  if(m.a->symm_=='G')      (*Vec) = m.a->eigenvector_r_;
  else if(m.a->symm_=='S') (*Vec) = m.a->eigenvector_;
  delete[] tmp;

  doTNT(this,Val,Vec);

  delete Vec; delete Val;
};
/* An old function, but I'd like to keep it around 
void Matrix<double>::operator=(const EXP& m){
  Matrix *Vec =  new Matrix(m.a->rows_,m.a->cols_);
  Matrix *Val =  new Matrix(m.a->rows_,m.a->cols_);
  double *tmp =  new double[m.a->rows_];
  double *tmp2 = new double[m.a->len_]; 
  Vec->clearAll();
  Val->clearAll();
//if(m.a->haveEigen_==0) { m.a->diag()};
  for(int i=0; i<m.a->rows_; i++){ tmp[i] = (double)exp(m.a->eigenvalue_[i]);};
  (*Vec) = (m.a->eigenvector_);
  Val->setDag(tmp);
  Vec->printAll();
  Val->printAll();
  double alpha = 1.0; double beta = 0.0;
  char xtrans = 'N'; char tmptrans = 'N';
  dgemm_(&xtrans,&Val->trans_,&Vec->rows_,&Val->cols_,&Vec->cols_,&alpha,Vec->data_,
         &Vec->rows_,Val->data_,&Vec->cols_,&beta,tmp2,&Vec->rows_);
  xtrans = 'T';
  dgemm_(&tmptrans,&xtrans,&Val->rows_,&Val->cols_,&Val->cols_,&alpha,tmp2,
	 &Val->rows_,Vec->data_,&Val->cols_,&beta,this->data_,&Val->rows_);
  delete[] tmp; delete[] tmp2;
  delete Vec; delete Val;
};
*/
} // namespace ChronusQ
