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

}; // QuasiNewton2<T>::buildSuperMatricies
