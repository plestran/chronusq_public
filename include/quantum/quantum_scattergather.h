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

template<>
template<typename Op>
void Quantum<double>::complexMyScale(Op &op){ };

template<>
template<typename Op>
void Quantum<dcomplex>::complexMyScale(Op &op){ op *= dcomplex(0.0,1.0); };

template<typename T>
void Quantum<T>::scatterDensity(){
  if(this->isScattered_) return;
  if(this->nTCS_ == 1 && this->isClosedShell)
    return;
  this->isScattered_ = true;

  // Allocate new scattered densities
  this->allocDensity(this->onePDMA_->cols() / this->nTCS_);

  std::vector<std::reference_wrapper<TMap>> scattered;
  scattered.emplace_back(*this->onePDMScalar_);
  scattered.emplace_back(*this->onePDMMz_);

  if(this->nTCS_ == 2) {
    scattered.emplace_back(*this->onePDMMy_);
    scattered.emplace_back(*this->onePDMMx_);
    Quantum<T>::spinScatter(*this->onePDMA_,scattered);
  } else if(!this->isClosedShell)
    Quantum<T>::spinScatter(*this->onePDMA_,*this->onePDMB_,scattered);

//this->onePDMA_.reset();
//this->onePDMB_.reset();
};

template<typename T>
void Quantum<T>::gatherDensity(){
  if(!this->isScattered_) return;
  if(this->nTCS_ == 1 && this->isClosedShell)
    return;
  this->isScattered_ = false;

  // Allocate new scattered densities
  this->allocDensity(this->onePDMScalar_->cols());

  std::vector<std::reference_wrapper<TMap>> scattered;
  scattered.emplace_back(*this->onePDMScalar_);
  scattered.emplace_back(*this->onePDMMz_);

  if(this->nTCS_ == 2) { 
    scattered.emplace_back(*this->onePDMMy_);
    scattered.emplace_back(*this->onePDMMx_);
    Quantum<T>::spinGather(*this->onePDMA_,scattered);
  } else if(!this->isClosedShell)
    Quantum<T>::spinGather(*this->onePDMA_,*this->onePDMB_,scattered);

  // Deallocate Space
//this->onePDMScalar_.reset();
//this->onePDMMz_.reset();
//this->onePDMMy_.reset();
//this->onePDMMx_.reset();
};

template<typename T>
template<typename Op>
void Quantum<T>::spinScatter(Op &op1, Op &op2, 
    std::vector<std::reference_wrapper<Op>> &scattered ){

  scattered[0].get() = op1 + op2;
  scattered[1].get() = op1 - op2;
}

template<typename T>
template<typename Op>
void Quantum<T>::spinScatter(Op &op, 
    std::vector<std::reference_wrapper<Op>> &scattered ){

  size_t currentDim = op.cols();

  Eigen::Map<TMatrix,0,Eigen::Stride<Dynamic,Dynamic> > 
    PAA(op.data(), currentDim/2, currentDim/2,
        Eigen::Stride<Dynamic,Dynamic>(2*currentDim,2));
  Eigen::Map<TMatrix,0,Eigen::Stride<Dynamic,Dynamic> > 
    PBA(op.data() + 1, currentDim/2, currentDim/2,
        Eigen::Stride<Dynamic,Dynamic>(2*currentDim,2));
  Eigen::Map<TMatrix,0,Eigen::Stride<Dynamic,Dynamic> > 
    PAB(op.data() + currentDim, currentDim/2, currentDim/2,
        Eigen::Stride<Dynamic,Dynamic>(2*currentDim,2));
  Eigen::Map<TMatrix,0,Eigen::Stride<Dynamic,Dynamic> > 
    PBB(op.data() + currentDim + 1, currentDim/2, currentDim/2,
        Eigen::Stride<Dynamic,Dynamic>(2*currentDim,2));

  scattered[0].get() = PAA + PBB;
  scattered[1].get() = PAA - PBB;
  scattered[2].get() = PAB - PBA;
  scattered[3].get() = PAB + PBA;

  // Scale My by "i"
  complexMyScale(scattered[2].get());

};

template<typename T>
template<typename Op>
void Quantum<T>::spinGather(Op &op1, Op &op2, 
    std::vector<std::reference_wrapper<Op>> &scattered ){

  op1.noalias() = 0.5 * (scattered[0].get() + scattered[1].get());
  op2.noalias() = 0.5 * (scattered[0].get() - scattered[1].get());
  
}

template<typename T>
template<typename Op>
void Quantum<T>::spinGather(Op &op, 
    std::vector<std::reference_wrapper<Op>> &scattered ){
  size_t currentDim = scattered[0].get().cols();

  // Since 
  //   My = i(PBA - PAB)
  //   PAB = Mx - i * My
  //   PBA = Mx + i * My
  // "My" is scaled by "i" before entering the reconstruction
  //
  // ** Note that this scaling is a dummy call for double 
  // precision objects and the sign is accounted for implicitly
  // through a flip in sign in the reconstruction **
  complexMyScale(scattered[2].get());

  Eigen::Map<TMatrix,0,Eigen::Stride<Dynamic,Dynamic> > 
    PAA(op.data(), currentDim, currentDim,
        Eigen::Stride<Dynamic,Dynamic>(2 * currentDim * 2,2));
  Eigen::Map<TMatrix,0,Eigen::Stride<Dynamic,Dynamic> > 
    PBA(op.data() + 1, currentDim, currentDim,
        Eigen::Stride<Dynamic,Dynamic>(2 * currentDim * 2,2));
  Eigen::Map<TMatrix,0,Eigen::Stride<Dynamic,Dynamic> > 
    PAB(op.data() + currentDim * 2,
        currentDim, currentDim,
        Eigen::Stride<Dynamic,Dynamic>(2 * currentDim * 2,2));
  Eigen::Map<TMatrix,0,Eigen::Stride<Dynamic,Dynamic> > 
    PBB(op.data() + (currentDim * 2) + 1,
        currentDim, currentDim,
        Eigen::Stride<Dynamic,Dynamic>(2 * currentDim * 2,2));

  PAA.noalias() = scattered[0].get() + scattered[1].get();
  PBB.noalias() = scattered[0].get() - scattered[1].get();

  if(typeid(T).hash_code() == typeid(dcomplex).hash_code()) {
    PAB.noalias() = scattered[3].get() - scattered[2].get();
    PBA.noalias() = scattered[3].get() + scattered[2].get();
  } else {
    // Sign flip viz complex case because there is an implied
    // "i" infront of the pure imaginary y component
    PAB.noalias() = scattered[3].get() + scattered[2].get();
    PBA.noalias() = scattered[3].get() - scattered[2].get();
  }
  op *= 0.5;
};

