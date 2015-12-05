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
void QuasiNewton2<T>::run(){
  time_t currentClockTime;

  time(&currentClockTime);

  this->out_ << "Quasi-Newton Calculation Started: " << 
                ctime(&currentClockTime) << endl;

  this->allocScr();
  auto start = std::chrono::high_resolution_clock::now();
  for(auto iter = 0; iter < this->maxMacroIter_; iter++){
    this->runMicro();
    if(this->isConverged_) break;
  };
  auto finish = std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> elapsed = finish - start;
  time(&currentClockTime);


}; // QuasiNewton2<T>::run

template<typename T>
void QuasiNewton2<T>::runMicro(){
  int NTrial = this->nGuess_;
  int NOld   = 0;
  int NNew   = this->nGuess_;


  this->readGuess();
  if(this->specialAlgorithm_ == SYMMETRIZED_TRIAL) this->symmetrizeTrial(); 

  for(auto iter = 0; iter < this->maxMicroIter_; iter++){
    this->formLinearTrans(NOld,NNew);
    CErr();
  };
}; // QuasiNewton2<T>::runMicro


template<typename T>
void QuasiNewton2<T>::formLinearTrans(const int NOld, const int NNew){
  auto N = this->qnObj_->nSingleDim();
  TMap NewSR  (this->SigmaRMem_ + (NOld*N),N,NNew);
  TMap NewVecR(this->TRMem_     + (NOld*N),N,NNew);

  TMap NewSL  (this->SigmaLMem_,0,0);
  TMap NewVecL(this->TLMem_    ,0,0);
  TMap NewRhoR(this->RhoRMem_  ,0,0);
  TMap NewRhoL(this->RhoLMem_  ,0,0);

  if(this->qnObj_->needsLeft()){
    new (&NewSL  ) TMap(this->SigmaLMem_ + (NOld*N),N,NNew);
    new (&NewVecL) TMap(this->TLMem_     + (NOld*N),N,NNew);
    new (&NewRhoR) TMap(this->RhoRMem_   + (NOld*N),N,NNew);
    new (&NewRhoL) TMap(this->RhoLMem_   + (NOld*N),N,NNew);
  }

  this->qnObj_->linearTrans(NewVecR,NewVecL,NewSR,NewSL,NewRhoR,NewRhoL);
};
