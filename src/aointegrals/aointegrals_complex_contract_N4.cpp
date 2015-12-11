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
  void AOIntegrals::twoEContractN4(bool RHF, bool KS, bool doFock, bool do24, bool doTCS,
    const ComplexMatrix &XAlpha, ComplexMatrix &AXAlpha, const ComplexMatrix &XBeta, 
    ComplexMatrix &AXBeta) {

    if(!this->haveAOTwoE) this->computeAOTwoE();
  
    if(getRank() == 0) {
      RealTensor2d ReXAlphaTensor, ReXBetaTensor;
      RealTensor2d ReAXAlphaTensor,ReAXBetaTensor;
      RealTensor2d ImXAlphaTensor, ImXBetaTensor;
      RealTensor2d ImAXAlphaTensor,ImAXBetaTensor;
      RealTensor2d ReXTotalTensor;
      RealTensor2d ImXTotalTensor;
     
      ReXAlphaTensor  = RealTensor2d(XAlpha.rows(),XAlpha.cols());
      ReAXAlphaTensor = RealTensor2d(AXAlpha.rows(),AXAlpha.cols());
      ImXAlphaTensor  = RealTensor2d(XAlpha.rows(),XAlpha.cols());
      ImAXAlphaTensor = RealTensor2d(AXAlpha.rows(),AXAlpha.cols());
      for(auto i = 0; i < XAlpha.size(); i++) {
        ReXAlphaTensor.storage()[i] = XAlpha.data()[i].real();
        ImXAlphaTensor.storage()[i] = XAlpha.data()[i].imag();
      }
      ReAXAlphaTensor.fill(0.0);
      ImAXAlphaTensor.fill(0.0);
      if(!RHF && !doTCS && !do24){
        ReXBetaTensor  = RealTensor2d(XBeta.rows(),XBeta.cols());
        ReAXBetaTensor = RealTensor2d(AXBeta.rows(),AXBeta.cols());
        ImXBetaTensor  = RealTensor2d(XBeta.rows(),XBeta.cols());
        ImAXBetaTensor = RealTensor2d(AXBeta.rows(),AXBeta.cols());
        for(auto i = 0; i < XBeta.size(); i++) {
          ReXBetaTensor.storage()[i] = XBeta.data()[i].real();
          ImXBetaTensor.storage()[i] = XBeta.data()[i].imag();
        }
        ReAXBetaTensor.fill(0.0);
        ImAXBetaTensor.fill(0.0);
     
        ReXTotalTensor = RealTensor2d(XAlpha.rows(),XBeta.cols());
        ImXTotalTensor = RealTensor2d(XAlpha.rows(),XBeta.cols());
        ReXTotalTensor = ReXAlphaTensor + ReXBetaTensor;
        ImXTotalTensor = ImXAlphaTensor + ImXBetaTensor;
      }
     
      enum{i,j,k,l}; 
   
      if(doFock)  {
        if(RHF){
          contract(1.0,*this->aoERI_,{i,j,k,l},ReXAlphaTensor,{l,k},0.0,
            ReAXAlphaTensor,{i,j});
          contract(1.0,*this->aoERI_,{i,j,k,l},ImXAlphaTensor,{l,k},0.0,
            ImAXAlphaTensor,{i,j});
          if(!KS){
            contract(-0.5,*this->aoERI_,{i,l,k,j},ReXAlphaTensor,{l,k},1.0,
              ReAXAlphaTensor,{i,j});
            contract(-0.5,*this->aoERI_,{i,l,k,j},ImXAlphaTensor,{l,k},1.0,
              ImAXAlphaTensor,{i,j});
          }
        } else if(!doTCS) {
          contract(1.0,*this->aoERI_,{i,j,k,l},ReXTotalTensor,{l,k},0.0,
            ReAXAlphaTensor,{i,j});
          contract(1.0,*this->aoERI_,{i,j,k,l},ImXTotalTensor,{l,k},0.0,
            ImAXAlphaTensor,{i,j});
          ReAXBetaTensor = ReAXAlphaTensor;
          ImAXBetaTensor = ImAXAlphaTensor;
          if(!KS) {
            contract(-1.0,*this->aoERI_,{i,l,k,j},ReXAlphaTensor,{l,k},1.0,
              ReAXAlphaTensor,{i,j});
            contract(-1.0,*this->aoERI_,{i,l,k,j},ImXAlphaTensor,{l,k},1.0,
              ImAXAlphaTensor,{i,j});
            contract(-1.0,*this->aoERI_,{i,l,k,j},ReXBetaTensor,{l,k},1.0,
              ReAXBetaTensor,{i,j});
            contract(-1.0,*this->aoERI_,{i,l,k,j},ImXBetaTensor,{l,k},1.0,
              ImAXBetaTensor,{i,j});
          }
        } else if(doTCS) {
          contract(1.0,*this->aoERI_,{i,j,k,l},ReXAlphaTensor,{l,k},0.0,
            ReAXAlphaTensor,{i,j});
          contract(1.0,*this->aoERI_,{i,j,k,l},ImXAlphaTensor,{l,k},0.0,
            ImAXAlphaTensor,{i,j});
          if(!KS) {
            contract(-1.0,*this->aoERI_,{i,l,k,j},ReXAlphaTensor,{l,k},1.0,
              ReAXAlphaTensor,{i,j});
            contract(-1.0,*this->aoERI_,{i,l,k,j},ImXAlphaTensor,{l,k},1.0,
              ImAXAlphaTensor,{i,j});
          }
        }
      } else if(do24) {
        contract(1.0,*this->aoERI_,{i,k,j,l},ReXAlphaTensor,{k,l},0.0,
          ReAXAlphaTensor,{i,j});
        contract(1.0,*this->aoERI_,{i,k,j,l},ImXAlphaTensor,{k,l},0.0,
          ImAXAlphaTensor,{i,j});
      }
      for(auto i = 0; i < AXAlpha.size(); i++) 
        AXAlpha.data()[i] = 
          dcomplex(ReAXAlphaTensor.storage()[i],ImAXAlphaTensor.storage()[i]);
      if(!RHF && !doTCS && !do24)
        for(auto i = 0; i < AXBeta.size(); i++) 
          AXBeta.data()[i] = 
            dcomplex(ReAXBetaTensor.storage()[i],ImAXBetaTensor.storage()[i]);
    } //serial code
#ifdef CQ_ENABLE_MPI
  // Syncronize MPI processes after PT build
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  }  // twoEContractN4
}; // namespace ChronusQ
