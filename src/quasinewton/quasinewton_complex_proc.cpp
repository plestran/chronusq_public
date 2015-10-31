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
  /** Reconstruct full dimenstion solution **/
  /*
   *  Reconstruct the approximate eigenvectors
   *
   *  For Hermitian matricies in general (viz. Davidson J. Comput. Phys. 17 (1975))
   *
   *  | X > = | b_i > X(R)_i
   *
   *  For RPA (viz. Kauczor et al. JCTC 7 (2010) p. 1610  (Eq 90,91)):
   *
   *  | X_g > = | {b_g}_i > {X(R)_g}_i
   *  | X_u > = | {b_u}_i > {X(R)_u}_i
   */ 
  template<>
  void QuasiNewton<dcomplex>::reconstructSolution(const int NTrial){
    ComplexCMMap XTSigmaR (this->XTSigmaRMem,NTrial,  NTrial);
    ComplexCMMap UR       (this->URMem,      this->N_,NTrial);
    ComplexCMMap TrialVecR(this->TVecRMem,   this->N_,NTrial);
    ComplexCMMap XTRhoR   (this->XTRhoRMem,  0,0);
    ComplexCMMap XTSigmaL (this->XTSigmaLMem,0,0);
    ComplexCMMap XTRhoL   (this->XTRhoLMem,  0,0);
    ComplexCMMap UL       (this->ULMem,      0,0);
    ComplexCMMap TrialVecL(this->TVecLMem,   0,0);
    RealVecMap ER(this->RealEMem,NTrial);
    if(!this->isHermitian_ || this->symmetrizedTrial_){
      new (&XTRhoR   ) ComplexCMMap(this->XTRhoRMem,  NTrial,  NTrial);
      new (&XTSigmaL ) ComplexCMMap(this->XTSigmaLMem,NTrial,  NTrial);
      new (&XTRhoL   ) ComplexCMMap(this->XTRhoLMem,  NTrial,  NTrial);
      new (&UL       ) ComplexCMMap(this->ULMem,      this->N_,NTrial);
      new (&TrialVecL) ComplexCMMap(this->TVecLMem,   this->N_,NTrial);
    }
    UR = TrialVecR * XTSigmaR;
    if(this->symmetrizedTrial_) UL = TrialVecL * XTSigmaL;
    // Stash away current approximation of eigenvalues and eigenvectors (NSek)
    (*this->solutionValues_) = ER.block(0,0,nSek_,1);
    (*this->solutionVector_) = UR.block(0,0,N_,nSek_); 
    // | X > = | X_g > + | X_u > (viz. Kauczor et al. JCTC 7 (2010) p. 1610  (Eq 80))
    if(this->symmetrizedTrial_)(*solutionVector_) += UL.block(0,0,N_,nSek_); 
  } // reconstructSolution

  /** Run Micro Iteration **/
  template<>
  void QuasiNewton<dcomplex>::runMicro(ostream &output){
    // Inital Values
    int NTrial = this->nGuess_;
    int NOld   = 0;
    int NNew   = this->nGuess_;
  
    // Initialize Trial Vectors
    ComplexCMMap TrialVecR(this->TVecRMem,this->N_,NTrial);
    ComplexCMMap TrialVecL(this->TVecLMem,0,0);
    if(!this->isHermitian_ || this->symmetrizedTrial_){
      new (&TrialVecL) ComplexCMMap(this->TVecLMem,this->N_,NTrial);
    }

    // Copy guess into Trial Vec
    TrialVecR = (*this->guessR_);
  //cout << TrialVecR << endl << endl;
    if(this->symmetrizedTrial_ || !this->isHermitian_){
      if(this->doRestart_) TrialVecL = (*this->guessL_);
      else                 TrialVecL = (*this->guessR_);
    }

    // Deallocate extraneous Guess storage
    this->guessR_.reset();
    if(this->doRestart_ && (this->symmetrizedTrial_ || !this->isHermitian_)) 
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
      this->nMicroIter_++;
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
