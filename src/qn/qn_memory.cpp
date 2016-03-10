/* 
 *  The Chronus Quantum (ChronusQ) software package is high-performace 
 *  computational chemistry software with a strong emphasis on explicitly 
 *  time-dependent and post-SCF quantum mechanical methods.
 *  
 *  Copyright (C) 2014-2016 Li Research Group (University of Washington)
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
#include <qn.h>

namespace ChronusQ{

template<>
void QuasiNewton2<double>::allocScrSpecial(){
  this->LRWORK = 0;
  if(this->problemType_ == DIAGONALIZATION){
    this->ERMem_ = new double[this->maxSubSpace_];
    if(this->matrixType_ != HERMETIAN)
      this->EIMem_ = new double[this->maxSubSpace_];
  }
};

template<>
void QuasiNewton2<dcomplex>::allocScrSpecial(){

  auto N      = this->qnObj_->nSingleDim();

  if(this->problemType_ == DIAGONALIZATION){
    if(this->matrixType_ == HERMETIAN)
      this->LRWORK = 2*N;
    else
      this->LWORK = 3*N;
    this->ECMem_ = new dcomplex[this->maxSubSpace_];
    this->ERMem_ = new   double[this->maxSubSpace_];
    this->RWORK_ = new   double[this->LRWORK];
  }
}

template<>
void QuasiNewton2<double>::cleanupScrSpecial(){
  if(this->problemType_ == DIAGONALIZATION){
    delete [] this->ERMem_;
    if(this->matrixType_ != HERMETIAN) delete [] this->EIMem_;
  }
};

template<>
void QuasiNewton2<dcomplex>::cleanupScrSpecial(){
  if(this->problemType_ == DIAGONALIZATION){
    delete [] this->ECMem_; 
    delete [] this->ERMem_; 
    delete [] this->RWORK_; 
  }
}
}; // namespace ChronusQ
