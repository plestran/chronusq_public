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
  /** Run Micro Iteration **/
  template<>
  void QuasiNewton<double>::runMicro(ostream &output){
    // Inital Values
    int NTrial = this->nGuess_;
    int NOld   = 0;
    int NNew   = this->nGuess_;
  
    // Initialize Trial Vectors
    RealCMMap TrialVecR(this->TVecRMem,this->N_,NTrial);
    RealCMMap TrialVecL(this->TVecLMem,0,0);
    if(!this->isHermetian_ || this->symmetrizedTrial_){
      new (&TrialVecL) RealCMMap(this->TVecLMem,this->N_,NTrial);
    }

    // Copy guess into Trial Vec
    TrialVecR = (*this->guessR_);
  //cout << TrialVecR << endl << endl;
    if(this->symmetrizedTrial_ || !this->isHermetian_){
      if(this->doRestart_) TrialVecL = (*this->guessL_);
      else                 TrialVecL = (*this->guessR_);
    }

    // Deallocate extraneous Guess storage
    this->guessR_.reset();
    if(this->doRestart_ && (this->symmetrizedTrial_ || !this->isHermetian_)) 
      this->guessL_.reset();

    // Symmetrize the trial vectors viz Kauczor et al. JCTC 7 (2010)
    if(!this->doRestart_ && this->symmetrizedTrial_){
      TrialVecR.block(this->N_/2,0,this->N_/2,this->nGuess_)
        = TrialVecR.block(0,0,this->N_/2,this->nGuess_);
      TrialVecL.block(this->N_/2,0,this->N_/2,this->nGuess_)
        = -TrialVecL.block(0,0,this->N_/2,this->nGuess_);
      // Normalize
      TrialVecR *= std::sqrt(0.5);
      TrialVecL *= std::sqrt(0.5);
    }

    for(auto iter = 0; iter < this->maxIter_; iter++){
      std::chrono::high_resolution_clock::time_point start,finish;
      std::chrono::duration<double> elapsed;
      output << "Starting Quasi-Newton Micro Iteration " << iter + 1 << endl;
      start = std::chrono::high_resolution_clock::now();
 
      // Perform the Linear Transformation
      this->linearTrans(NOld,NNew);
      // Perform full projection
      this->fullProjection(NTrial);
      // Build Supermatricies if symmetrized vectors are used
      if(this->symmetrizedTrial_) this->buildSuperMat(NTrial);
      // Diagonalize the subspace
      if(this->doDiag_) this->redDiag(NTrial,output);
      // Reconstruct the solution vectors
      this->reconstructSolution(NTrial); 
      // Generate the residual vectors
      this->genRes(NTrial);
      // Check the convergence
      int NNotConv;
      auto resConv = this->checkConv(NTrial,NNotConv,output); 
      // If NSek (lowest) residuals are converged, exit, else, create
      // perturbed guess vectors for next iteration
      this->isConverged_ = (NNotConv == 0);
      this->doRestart_   = (NTrial+NNotConv > this->maxSubSpace_);
      if(this->isConverged_) {
        finish = std::chrono::high_resolution_clock::now();
        elapsed = finish - start;
        output << "  Quasi-Newton Micro Iteration took " << std::fixed 
               << elapsed.count() << " secs" << endl << endl;
        break;
      }
      if(this->doRestart_) {
        this->setupRestart();
        break;
      }
      this->formNewGuess(resConv,NTrial,NNotConv,NOld,NNew);

      finish = std::chrono::high_resolution_clock::now();
      elapsed = finish - start;
      output << "Quasi-Newton Micro Iteration took " << std::fixed 
             << elapsed.count() << " secs" << endl << endl;

    } // for iter in [0,maxIter)
  }; // runMicro
}; // namespace ChronusQ
