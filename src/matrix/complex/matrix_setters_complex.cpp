#include "matrix.h"
using ChronusQ::Matrix;

namespace ChronusQ {
template<>
template<> void Matrix<dcomplex>::cpyData<dcomplex>(dcomplex *x){ 
  for(int i = 0; i < this->len_; i++){
    this->data_[i] = x[i];
  };
};
template<>
template<> void Matrix<dcomplex>::setDag<dcomplex>(dcomplex *x){
  if(this->format_!=0){ throw 3020;}
  else {
    for(int i=0;i<this->cols_;i++){ this->data_[this->rows_*i + i] = x[i];};
  };
}
template<>
template<> void Matrix<dcomplex>::setDag<double>(double *x){
  if(this->format_!=0){ throw 3020;}
  else {
    for(int i=0;i<this->cols_;i++){ this->data_[this->rows_*i + i] = x[i];};
  };
}
} // namespace ChronusQ
