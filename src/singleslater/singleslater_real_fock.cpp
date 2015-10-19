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
#include <singleslater.h>
namespace ChronusQ {
#ifndef USE_LIBINT
//----------------------------//
// form the Coulomb matrix    //
//----------------------------//
template<>
void SingleSlater<double>::formCoulomb(){
//clock_t start,finish;
//start = clock();
  std::chrono::high_resolution_clock::time_point start, finish;
  start = std::chrono::high_resolution_clock::now();
  int i;
  
  if(!this->haveDensity) this->formDensity();
  if(!this->aointegrals_->haveAOTwoE) this->aointegrals_->computeAOTwoE();
  this->coulombA_->setZero();

  for(i=0;i<this->nBasis_;i++) (*densityA_)(i,i)*=math.half;

  this->densityA_->vectorize();
  this->coulombA_->vectorize();
  (*this->coulombA_) = (*this->aointegrals_->twoEC_)*(*this->densityA_);
  this->coulombA_->unvectorize();
  this->densityA_->unvectorize();

  for(i=0;i<this->nBasis_;i++) (*densityA_)(i,i)*=math.two;

  finish = std::chrono::high_resolution_clock::now();
  this->aointegrals_->CoulD = finish - start; 
  this->fileio_->out<<"\nCPU time for building the Coulomb matrix:  "<< this->aointegrals_->CoulD.count() <<" seconds."<<endl;

//if(this->printLevel_ >= 2) this->coulombA_->printAll(5,this->fileio_->out);
  if(this->printLevel_ >= 2) prettyPrint(this->fileio_->out,(*this->coulombA_),"Alpha Coulomb");
  this->haveCoulomb = true;
};
//----------------------------//
// form the exchange matrix    //
//----------------------------//
template<>
void SingleSlater<double>::formExchange(){
//clock_t start,finish;
//start = clock();
  std::chrono::high_resolution_clock::time_point start, finish;
  start = std::chrono::high_resolution_clock::now();
  int i;
  
  if(!this->haveDensity) this->formDensity();
  if(!this->aointegrals_->haveAOTwoE) this->aointegrals_->computeAOTwoE();
  this->exchangeA_->setZero();

  for(i=0;i<this->nBasis_;i++) (*densityA_)(i,i)*=math.half;

  this->densityA_->vectorize();
  this->exchangeA_->vectorize();
  (*this->exchangeA_)= (*this->aointegrals_->twoEX_)*(*this->densityA_);
  this->exchangeA_->scale(math.quarter);
  this->exchangeA_->unvectorize();
  this->densityA_->unvectorize();

  for(i=0;i<this->nBasis_;i++) (*densityA_)(i,i)*=math.two;

  finish = std::chrono::high_resolution_clock::now();
  this->aointegrals_->ExchD = finish - start; 
  this->fileio_->out<<"\nCPU time for building the Exchange matrix:  "<< this->aointegrals_->ExchD.count() <<" seconds."<<endl;

//if(this->printLevel_ >= 2) this->exchangeA_->printAll(5,this->fileio_->out);
  if(this->printLevel_ >= 2) prettyPrint(this->fileio_->out,(*this->exchangeA_),"Alpha Exchange");

  this->haveExchange = true;
};
#endif
}; // namespace ChronusQ
