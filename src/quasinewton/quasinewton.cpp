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
void QuasiNewton<double>::redDiag(int NTrial,ostream &output){
  // LAPACK Variables
  char JOBVR = 'V'; // Get Right eigenvectors
  char JOBVL = 'N'; // Get Left eigenvectors
  char UPLO = 'L';  // Lower triagnle is used for symmetrix matricies
  int INFO;         // Success / failure flag
  int iType = 1;    // Flag for the type of GEP being solved
  int TwoNTrial = 2*NTrial; // Dimension of the supermatricies

  int IOff = NTrial; // Space for eigenvalues
  if(!this->isHermetian_ || this->symmetrizedTrial_){
    IOff += NTrial; // Space for paired eigenvalues or imaginary part
    new (&this->ER) RealVecMap(this->LAPACK_SCR,2*NTrial);
  } else{ new (&this->ER) RealVecMap(this->LAPACK_SCR,NTrial);}
  double * WORK = this->LAPACK_SCR + IOff;

  if(this->isHermetian_) {
    if(!this->symmetrizedTrial_){
      // Solve E(R)| X(R) > = | X(R) > ω
      dsyev_(&JOBVR,&UPLO,&NTrial,this->XTSigmaR.data(),&NTrial,
             this->ER.data(),WORK,&this->LWORK,&INFO); 
      if(INFO!=0) CErr("DSYEV failed to converge in Davison Iterations",output);
    } else {
      RealCMMatrix SCPY(SSuper); // Copy of original matrix to use for re-orthogonalization
      /*
       * Solve S(R)| X(R) > = E(R)| X(R) > (1/ω)
       *
       * | X(R) > = | X(R)_g >
       *            | X(R)_u >
       *
       * The opposite (1/ω vs ω) is solved because the metric is not positive definite
       * and can therefore not be solved using DSYGV because of the involved Cholesky
       * decomposition.
       *
       * This must be revisited for Spinor orbital hessians! FIXME
       */ 
      dsygv_(&iType,&JOBVR,&UPLO,&TwoNTrial,this->SSuper.data(),&TwoNTrial,
             this->ASuper.data(),&TwoNTrial,this->ER.data(),WORK,&this->LWORK,
             &INFO);
      if(INFO!=0) CErr("DSYGV failed to converge in Davison Iterations",output);

      // Grab the "positive paired" roots (throw away other element of the pair)
      new (&this->ER)     RealVecMap(this->LAPACK_SCR+NTrial,NTrial);
      new (&this->SSuper) RealCMMap (this->SSuperMem+2*NTrial*NTrial,2*NTrial,NTrial);

      // Swap the ordering because we solve for (1/ω)
      for(auto i = 0 ; i < NTrial; i++) ER(i) = 1.0/this->ER(i);
      for(auto i = 0 ; i < NTrial/2; i++){
        this->SSuper.col(i).swap(this->SSuper.col(NTrial - i - 1));
        double tmp = this->ER(i);
        this->ER(i) = this->ER(NTrial - i - 1);
        this->ER(NTrial - i - 1) = tmp;
      }

      // Re-orthogonalize the eigenvectors with respect to the metric S(R)
      // because DSYGV orthogonalzies the vectors with respect to E(R)
      // because we solve the opposite problem.
      //
      // Gramm-Schmidt
      for(auto i = 0; i < NTrial; i++){
        double inner = this->SSuper.col(i).dot(SCPY*this->SSuper.col(i));
        int sgn = inner / std::abs(inner);
        inner = sgn*std::sqrt(sgn*inner);
        this->SSuper.col(i) /= inner;
        for(auto j = i+1; j < NTrial; j++){
          this->SSuper.col(j) -= 
            this->SSuper.col(i)*
            (this->SSuper.col(i).dot(SCPY*this->SSuper.col(j)));
        }
      }

      this->XTSigmaR = this->SSuper.block(0,     0,NTrial,NTrial);
      this->XTSigmaL = this->SSuper.block(NTrial,0,NTrial,NTrial);
    }
  } 
} // redDiag
template<>
void QuasiNewton<double>::runMicro(ostream &output){
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
    // Perform the Linear Transformation
    this->linearTrans();
    // Perform full projection
    this->fullProjection();
    // Build Supermatricies if symmetrized vectors are used
    if(this->symmetrizedTrial_) this->buildSuperMat(NTrial);
    // Diagonalize the subspace
    if(this->doDiag_) this->redDiag(NTrial);
    // Reconstruct the solution vectors
    this->reconstructSolution(); 

  } // for iter in [0, maxIter)
  
} // runMicro

} // namespace ChronusQ
