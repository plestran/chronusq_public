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
