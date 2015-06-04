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
#include <aointegrals.h>
namespace ChronusQ{
#ifdef USE_LIBINT
template<>
void AOIntegrals::twoEContractDirect(bool RHF, bool doFock, const ComplexMatrix &XAlpha, ComplexMatrix &AXAlpha,
                                 const ComplexMatrix &XBeta, ComplexMatrix &AXBeta) {
  CErr("No Direct Contraction for Complex Matricies Implemented");
}
template<>
void AOIntegrals::twoEContractN4(bool RHF, bool doFock, const ComplexMatrix &XAlpha, ComplexMatrix &AXAlpha,
                                 const ComplexMatrix &XBeta, ComplexMatrix &AXBeta) {
  this->fileio_->out << "Contracting with in-core two-electron integrals" << endl;
  if(!this->haveAOTwoE) this->computeAOTwoE();
  RealTensor2d ReXAlphaTensor(XAlpha.rows(),XAlpha.cols());
  RealTensor2d ImXAlphaTensor(XAlpha.rows(),XAlpha.cols());
  RealTensor2d ReAXAlphaTensor(AXAlpha.rows(),AXAlpha.cols());
  RealTensor2d ImAXAlphaTensor(AXAlpha.rows(),AXAlpha.cols());
  for(auto i = 0; i < XAlpha.size(); i++) {
    ReXAlphaTensor.storage()[i] = XAlpha.data()[i].real();
    ImXAlphaTensor.storage()[i] = XAlpha.data()[i].imag();
  }
  ReAXAlphaTensor.fill(0.0);
  ImAXAlphaTensor.fill(0.0);

  double fact = -1.0;
  if(RHF && doFock) fact = -0.5;

  enum{i,j,k,l}; 
  contract(1.0,*this->aoERI_,{i,j,k,l},ReXAlphaTensor,{l,k},0.0,ReAXAlphaTensor,{i,j});
  contract(1.0,*this->aoERI_,{i,j,k,l},ImXAlphaTensor,{l,k},0.0,ImAXAlphaTensor,{i,j});
  contract(fact,*this->aoERI_,{i,l,k,j},ReXAlphaTensor,{l,k},1.0,ReAXAlphaTensor,{i,j});
  contract(fact,*this->aoERI_,{i,l,k,j},ImXAlphaTensor,{l,k},1.0,ImAXAlphaTensor,{i,j});

  for(auto i = 0; i < AXAlpha.size(); i++) {
    AXAlpha.data()[i] = 
      dcomplex(ReAXAlphaTensor.storage()[i],ImAXAlphaTensor.storage()[i]);
  }
}
template<>
void AOIntegrals::twoEContractDF(bool RHF, bool doFock, const ComplexMatrix &XAlpha, ComplexMatrix &AXAlpha,
                                 const ComplexMatrix &XBeta, ComplexMatrix &AXBeta) {
  CErr("No Density Fitting Contraction for Complex Matricies Implemented");
}
#endif

}
