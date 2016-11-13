/* 
 *  This file is part of the Chronus Quantum (ChronusQ) software package
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
template<typename Op,
  typename std::enable_if<
    std::is_same<typename Op::Scalar,
                 double>::value,int>::type>
void AOIntegrals::Ortho1Trans(Op& op1, Op& op2){
  double* Scratch = 
    this->memManager_->malloc<double>(op1.size());
  RealMap SCR(Scratch,op1.rows(),op1.rows());

  SCR.noalias() = (*this->ortho1_) * op1;
  op2.noalias() = SCR * this->ortho1_->transpose();
  

  this->memManager_->free(Scratch,op1.size());
}

template<typename Op,
  typename std::enable_if<
    std::is_same<typename Op::Scalar,
                 dcomplex>::value,int>::type>
void AOIntegrals::Ortho1Trans(Op& op1, Op& op2){
  dcomplex* Scratch = 
    this->memManager_->malloc<dcomplex>(op1.size());
  ComplexMap SCR(Scratch,op1.rows(),op1.rows());

  SCR.real() = (*this->ortho1_) * op1.real();
  SCR.imag() = (*this->ortho1_) * op1.imag();
  op2.real() = SCR.real() * this->ortho1_->transpose();
  op2.imag() = SCR.imag() * this->ortho1_->transpose();
  

  this->memManager_->free(Scratch,op1.size());
}

template<typename Op,
  typename std::enable_if<
    std::is_same<typename Op::Scalar,
                 double>::value,int>::type>
void AOIntegrals::Ortho1TransT(Op& op1, Op& op2){
  double* Scratch = 
    this->memManager_->malloc<double>(op1.size());
  RealMap SCR(Scratch,op1.rows(),op1.rows());

  SCR.noalias() = this->ortho1_->transpose() * op1;
  op2.noalias() = SCR * (*this->ortho1_);
  

  this->memManager_->free(Scratch,op1.size());
}

template<typename Op,
  typename std::enable_if<
    std::is_same<typename Op::Scalar,
                 dcomplex>::value,int>::type>
void AOIntegrals::Ortho1TransT(Op& op1, Op& op2){
  dcomplex* Scratch = 
    this->memManager_->malloc<dcomplex>(op1.size());
  ComplexMap SCR(Scratch,op1.rows(),op1.rows());

  SCR.real() = this->ortho1_->transpose() * op1.real();
  SCR.imag() = this->ortho1_->transpose() * op1.imag();
  op2.real() = SCR.real() * (*this->ortho1_);
  op2.imag() = SCR.imag() * (*this->ortho1_);
  

  this->memManager_->free(Scratch,op1.size());
}


template<typename Op,
  typename std::enable_if<
    std::is_same<typename Op::Scalar,
                 double>::value,int>::type>
void AOIntegrals::Ortho2Trans(Op& op1, Op& op2){
  double* Scratch = 
    this->memManager_->malloc<double>(op1.size());
  RealMap SCR(Scratch,op1.rows(),op1.rows());

  SCR.noalias() = this->ortho2_->transpose() * op1;
  op2.noalias() = SCR * (*this->ortho2_);
  

  this->memManager_->free(Scratch,op1.size());
}

template<typename Op,
  typename std::enable_if<
    std::is_same<typename Op::Scalar,
                 dcomplex>::value,int>::type>
void AOIntegrals::Ortho2Trans(Op& op1, Op& op2){
  dcomplex* Scratch = 
    this->memManager_->malloc<dcomplex>(op1.size());
  ComplexMap SCR(Scratch,op1.rows(),op1.rows());

  SCR.real() = this->ortho2_->transpose() * op1.real();
  SCR.imag() = this->ortho2_->transpose() * op1.imag();
  op2.real() = SCR.real() * (*this->ortho2_);
  op2.imag() = SCR.imag() * (*this->ortho2_);
  

  this->memManager_->free(Scratch,op1.size());
}

