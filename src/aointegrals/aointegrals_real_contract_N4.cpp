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
  void AOIntegrals::twoEContractN4(bool RHF, bool doFock, bool doTCS,
    const RealMatrix &XAlpha, RealMatrix &AXAlpha, const RealMatrix &XBeta, 
    RealMatrix &AXBeta) {

    this->fileio_->out << "Contracting with in-core two-electron integrals" << endl;
    if(!this->haveAOTwoE) this->computeAOTwoE();
  
    RealTensor2d XAlphaTensor,XBetaTensor;
    RealTensor2d AXAlphaTensor,AXBetaTensor;
    enum{i,j,k,l}; 
  
    if(doFock)  {
      XAlphaTensor  = RealTensor2d(XAlpha.rows(),XAlpha.cols());
      AXAlphaTensor = RealTensor2d(AXAlpha.rows(),AXAlpha.cols());
      for(auto i = 0; i < XAlpha.size(); i++) XAlphaTensor.storage()[i] = XAlpha.data()[i];
      AXAlphaTensor.fill(0.0);
      if(RHF){
        contract(1.0,*this->aoERI_,{i,j,k,l},XAlphaTensor,{l,k},0.0,AXAlphaTensor,{i,j});
        contract(-0.5,*this->aoERI_,{i,l,k,j},XAlphaTensor,{l,k},1.0,AXAlphaTensor,{i,j});
      } else if(!doTCS) {
        XBetaTensor  = RealTensor2d(XBeta.rows(),XBeta.cols());
        AXBetaTensor = RealTensor2d(AXBeta.rows(),AXBeta.cols());
        for(auto i = 0; i < XBeta.size(); i++) XBetaTensor.storage()[i] = XBeta.data()[i];
        AXBetaTensor.fill(0.0);
        RealTensor2d XTotalTensor(XAlpha.rows(),XBeta.cols());
  
        XTotalTensor = XAlphaTensor + XBetaTensor;
  
        contract(1.0,*this->aoERI_,{i,j,k,l},XTotalTensor,{l,k},0.0,AXAlphaTensor,{i,j});
        AXBetaTensor = AXAlphaTensor;
        contract(-1.0,*this->aoERI_,{i,l,k,j},XAlphaTensor,{l,k},1.0,AXAlphaTensor,{i,j});
        contract(-1.0,*this->aoERI_,{i,l,k,j},XBetaTensor,{l,k},1.0,AXBetaTensor,{i,j});
      } else if(doTCS) {
        contract(1.0 ,*this->aoERI_,{i,j,k,l},XAlphaTensor,{l,k},0.0,AXAlphaTensor,{i,j});
        contract(-1.0,*this->aoERI_,{i,l,k,j},XAlphaTensor,{l,k},1.0,AXAlphaTensor,{i,j});
      }
    } else CErr("General Contraction NYI for in-core integrals");
     for(auto i = 0; i < AXAlpha.size(); i++) AXAlpha.data()[i] = AXAlphaTensor.storage()[i];
     if(!RHF && !doTCS)
       for(auto i = 0; i < AXBeta.size(); i++) AXBeta.data()[i] = AXBetaTensor.storage()[i];
  }  // twoEContractN4
}; // namespace ChronusQ
