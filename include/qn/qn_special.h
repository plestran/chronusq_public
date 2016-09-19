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

template<typename T>
void QuasiNewton2<T>::symmetrizeTrial(){
  auto N      = this->qnObj_->nSingleDim();
  auto NGuess = this->qnObj_->nGuess();
  TMap TVecR   (this->TRMem_, N, NGuess);
  TMap TVecL   (this->TLMem_, N, NGuess);

  TVecR.block(N/2,0,N/2,NGuess) =  TVecR.block(0,0,N/2,NGuess);
  TVecL.block(N/2,0,N/2,NGuess) = -TVecL.block(0,0,N/2,NGuess);

  TVecR *= std::sqrt(0.5);
  TVecL *= std::sqrt(0.5);


}; // QuasiNewton2<T>::symmetrizeTrial

template<typename T>
void QuasiNewton2<T>::buildSuperMatricies(const int NTrial){
  TMap XTSigmaR(this->XTSigmaRMem_,NTrial,  NTrial);
  TMap XTRhoR  (this->XTRhoRMem_,  NTrial,  NTrial);
  TMap XTSigmaL(this->XTSigmaLMem_,NTrial,  NTrial);
  TMap XTRhoL  (this->XTRhoLMem_,  NTrial,  NTrial);
  TMap ASuper  (this->ASuperMem_, 2*NTrial,2*NTrial);
  TMap SSuper  (this->SSuperMem_, 2*NTrial,2*NTrial);

  ASuper.setZero();
  SSuper.setZero();
  ASuper.block(0,     0,     NTrial,NTrial) = XTSigmaR;
  ASuper.block(NTrial,NTrial,NTrial,NTrial) = XTSigmaL;
  SSuper.block(0,     NTrial,NTrial,NTrial) = XTRhoR;
  SSuper.block(NTrial,0,     NTrial,NTrial) = XTRhoL;

  hsize_t offset[] = {0,0};
  hsize_t stride[] = {1,1};
  hsize_t block[]  = {1,1};

  hsize_t subDim[] = {2*NTrial,2*NTrial};
  hsize_t count[]  = {2*NTrial,2*NTrial};
 
  H5::DataSpace memSpace(2,subDim,NULL);
  H5::DataSpace subDataSpace;
  subDataSpace = this->ASuperFile_->getSpace();
  subDataSpace.selectHyperslab(H5S_SELECT_SET,count,offset,stride,block);

  this->ASuperFile_->write(this->ASuperMem_,H5PredType<T>(),memSpace,
    subDataSpace);
  this->SSuperFile_->write(this->SSuperMem_,H5PredType<T>(),memSpace,
    subDataSpace);
//prettyPrint(cout,SSuper*SSuper,"SS");
//prettyPrint(cout,SSuper.inverse()*SSuper,"SinvS");
}; // QuasiNewton2<T>::buildSuperMatricies

template<typename T>
void QuasiNewton2<T>::formNHrProd(const int NTrial) {
  auto TwoNTrial = 2 * NTrial;
  this->invertSuperMetric(NTrial);


  TMap  SSuper(this->SSuperMem_, TwoNTrial,TwoNTrial);
  TMap  ASuper(this->ASuperMem_, TwoNTrial,TwoNTrial);
  TMap NHrProd(this->NHrProdMem_,TwoNTrial,TwoNTrial);

  NHrProd = SSuper * ASuper;
//NHrProd = SSuper.inverse() * ASuper;
}; // QuasiNewton2<T>::formNHrProd
