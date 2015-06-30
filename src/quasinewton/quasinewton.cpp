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
#include <quasinewton.h>

namespace ChronusQ {
template<>
void QuasiNewton<double>::runMicro(){
  // LAPACK Variables
  char JOBVR = 'V';
  char JOBVL = 'N';
  char UPLO = 'L';
  int INFO;

  // Inital Values
  int NTrial = this->nGuess_;
  int NOld   = 0;
  int NNew   = this->nGuess_;

  // Initialize Trial Vectors
  new (&this->TrialVecR) RealCMMap(this->TVecRMem,this->N_,NTrial);
  if(!this->isHermetian_ || this->symmetrizedTrial_){
    new (&this->TrialVecL) RealCMMap(this->TVecLMem,this->N_,NTrial);
  }

  // Symmetrize the trial vectors viz Kauczor et al. JCTC 7 (2010)
  if(this->symmetrizedTrial_){
    this->TrialVecR = (*this->guess_);
    this->TrialVecL = (*this->guess_);
    this->TrialVecR.block(this->N_/2,0,this->N_/2,this->nGuess_)
      = this->TrialVecR.block(0,0,this->N_/2,this->nGuess_);
    this->TrialVecL.block(this->N_/2,0,this->N_/2,this->nGuess_)
      = -this->TrialVecL.block(0,0,this->N_/2,this->nGuess_);
    // Normalize
    this->TrialVecR *= std::sqrt(0.5);
    this->TrialVecL *= std::sqrt(0.5);

  } else {
    this->TrialVecR = (*this->guess_);
  }

  for(auto iter = 0; iter < this->maxIter_; iter++){
    std::chrono::high_resolution_clock::time_point start,finish;
    std::chrono::duration<double> elapsed;
    output << "Starting Quasi-Newton Micro Iteration " << iter + 1 << endl;
    start = std::chrono::high_resolution_clock::now();

    // Resize the Eigen Mapes to fit new vectors
    this->resizeMaps(NTrial,NOld,NNew);

  } // for iter in [0, maxIter)
  
} // runMicro
} // namespace ChronusQ
