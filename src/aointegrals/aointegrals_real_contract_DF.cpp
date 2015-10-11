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
  template<>
  void AOIntegrals::twoEContractDF(bool RHF, bool KS, bool doFock, const RealMatrix &XAlpha, 
    RealMatrix &AXAlpha, const RealMatrix &XBeta, RealMatrix &AXBeta) {

//  this->fileio_->out << "Contracting with in-core density fitting integrals" << endl;
    if(!this->haveRIS)  this->computeAORIS();
    if(!this->haveRII)  this->computeAORII();
    if(!this->haveTRII) this->transformAORII();
  
    RealTensor2d XAlphaTensor(XAlpha.rows(),XAlpha.cols());
    RealTensor2d AXAlphaTensor(AXAlpha.rows(),AXAlpha.cols());
  
    for(auto i = 0; i < XAlpha.size(); i++) XAlphaTensor.storage()[i] = XAlpha.data()[i];
    AXAlphaTensor.fill(0.0);
  
    RealTensor1d A(this->DFbasisSet_->nBasis());
    RealTensor3d B(this->nBasis_,this->nBasis_,this->DFbasisSet_->nBasis());
  
    double fact = -1.0;
    if(RHF && doFock) fact = -0.5;
  
    enum{i,j,k,l,XX,YY}; 
    contract(1.0,*this->aoRII_,{i,j,XX},XAlphaTensor,{i,j},0.0,A,{XX});
    contract(1.0,*this->aoRII_,{k,i,XX},XAlphaTensor,{k,j},0.0,B,{i,j,XX});
    contract(1.0,*this->aoRII_,{i,j,XX},A,{XX},0.0,AXAlphaTensor,{i,j});
    contract(-0.5,*this->aoRII_,{i,k,XX},B,{j,k,XX},1.0,AXAlphaTensor,{i,j});
  
    for(auto i = 0; i < AXAlpha.size(); i++) AXAlpha.data()[i] = AXAlphaTensor.storage()[i];
  } // twoEContractDF

} // namespace ChronusQ

