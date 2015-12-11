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
    const RealMatrix &XAlpha, RealMatrix &AXAlpha, const RealMatrix &XBeta, 
    RealMatrix &AXBeta) {
#ifdef CQ_ENABLE_MPI
  // Syncronize MPI processes after PT build
  MPI_Barrier(MPI_COMM_WORLD);
  printf("Hello 5 %d:%d\n",getRank(),getSize());
#endif


    if(!this->haveAOTwoE) this->computeAOTwoE();
    if(getRank() == 0) {
   
      RealTensor2d XAlphaTensor,XBetaTensor;
      RealTensor2d AXAlphaTensor,AXBetaTensor;
      RealTensor2d XTotalTensor;
     
      XAlphaTensor  = RealTensor2d(XAlpha.rows(),XAlpha.cols());
      AXAlphaTensor = RealTensor2d(AXAlpha.rows(),AXAlpha.cols());
      for(auto i = 0; i < XAlpha.size(); i++) 
        XAlphaTensor.storage()[i] = XAlpha.data()[i];
      AXAlphaTensor.fill(0.0);
      if(!RHF && !doTCS && !do24){
        XBetaTensor  = RealTensor2d(XBeta.rows(),XBeta.cols());
        AXBetaTensor = RealTensor2d(AXBeta.rows(),AXBeta.cols());
        for(auto i = 0; i < XBeta.size(); i++) 
          XBetaTensor.storage()[i] = XBeta.data()[i];
        AXBetaTensor.fill(0.0);
     
        XTotalTensor = RealTensor2d(XAlpha.rows(),XBeta.cols());
        XTotalTensor = XAlphaTensor + XBetaTensor;
      }
     
      enum{i,j,k,l}; 
   
      if(doFock)  {
        if(RHF){
          contract(1.0,*this->aoERI_,{i,j,k,l},XAlphaTensor,{l,k},0.0,
            AXAlphaTensor,{i,j});
          if(!KS) contract(-0.5,*this->aoERI_,{i,l,k,j},XAlphaTensor,{l,k},1.0,
                    AXAlphaTensor,{i,j});
        } else if(!doTCS) {
          contract(1.0,*this->aoERI_,{i,j,k,l},XTotalTensor,{l,k},0.0,
            AXAlphaTensor,{i,j});
          AXBetaTensor = AXAlphaTensor;
          if(!KS) {
            contract(-1.0,*this->aoERI_,{i,l,k,j},XAlphaTensor,{l,k},1.0,
              AXAlphaTensor,{i,j});
            contract(-1.0,*this->aoERI_,{i,l,k,j},XBetaTensor,{l,k},1.0,
              AXBetaTensor,{i,j});
          }
        } else if(doTCS) {
          contract(1.0 ,*this->aoERI_,{i,j,k,l},XAlphaTensor,{l,k},0.0,
            AXAlphaTensor,{i,j});
          if(!KS) contract(-1.0,*this->aoERI_,{i,l,k,j},XAlphaTensor,{l,k},1.0,
                    AXAlphaTensor,{i,j});
        }
      } else if(do24) {
        contract(1.0,*this->aoERI_,{i,k,j,l},XAlphaTensor,{k,l},0.0,
          AXAlphaTensor,{i,j});
      }
      for(auto i = 0; i < AXAlpha.size(); i++) 
        AXAlpha.data()[i] = AXAlphaTensor.storage()[i];
      if(!RHF && !doTCS && !do24)
        for(auto i = 0; i < AXBeta.size(); i++) 
          AXBeta.data()[i] = AXBetaTensor.storage()[i];
    } // serial code
#ifdef CQ_ENABLE_MPI
  // Syncronize MPI processes after PT build
  MPI_Barrier(MPI_COMM_WORLD);
  printf("Hello 1 %d:%d\n",getRank(),getSize());
#endif
  }  // twoEContractN4
}; // namespace ChronusQ
