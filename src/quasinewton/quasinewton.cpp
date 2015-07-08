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
  /** Diagoanlize Reduced Problem **/
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
    }
    this->ERMem = LAPACK_SCR;
    this->WORK  = this->ERMem + IOff;
    if(this->isHermetian_) {
      if(!this->symmetrizedTrial_){
        // Solve E(R)| X(R) > = | X(R) > ω
        dsyev_(&JOBVR,&UPLO,&NTrial,this->XTSigmaRMem,&NTrial,
               this->ERMem,this->WORK,&this->LWORK,&INFO); 
        if(INFO!=0) CErr("DSYEV failed to converge in Davison Iterations",output);
      } else {
        RealCMMap SSuper(this->SSuperMem, 2*NTrial,2*NTrial);
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
        dsygv_(&iType,&JOBVR,&UPLO,&TwoNTrial,this->SSuperMem,&TwoNTrial,
               this->ASuperMem,&TwoNTrial,this->ERMem,this->WORK,&this->LWORK,
               &INFO);
        if(INFO!=0) CErr("DSYGV failed to converge in Davison Iterations",output);
        // Grab the "positive paired" roots (throw away other element of the pair)
        this->ERMem += NTrial;
        RealVecMap ER    (this->ERMem,NTrial);
        new (&SSuper) RealCMMap(this->SSuperMem+2*NTrial*NTrial,2*NTrial,NTrial);
 
        // Swap the ordering because we solve for (1/ω)
        for(auto i = 0 ; i < NTrial; i++) ER(i) = 1.0/ER(i);
        for(auto i = 0 ; i < NTrial/2; i++){
          SSuper.col(i).swap(SSuper.col(NTrial - i - 1));
          double tmp = ER(i);
          ER(i) = ER(NTrial - i - 1);
          ER(NTrial - i - 1) = tmp;
        }
 
        // Re-orthogonalize the eigenvectors with respect to the metric S(R)
        // because DSYGV orthogonalzies the vectors with respect to E(R)
        // because we solve the opposite problem.
        //
        // Gramm-Schmidt
        for(auto i = 0; i < NTrial; i++){
          double inner = SSuper.col(i).dot(SCPY*SSuper.col(i));
          int sgn = inner / std::abs(inner);
          inner = sgn*std::sqrt(sgn*inner);
          SSuper.col(i) /= inner;
          for(auto j = i+1; j < NTrial; j++){
            SSuper.col(j) -= 
              SSuper.col(i)*(SSuper.col(i).dot(SCPY*SSuper.col(j)));
          }
        }
 
        // Separate the eigenvectors into gerade and ungerade parts
        RealCMMap XTSigmaR(this->XTSigmaRMem,NTrial,NTrial);
        RealCMMap XTSigmaL(this->XTSigmaLMem,NTrial,NTrial);
        XTSigmaR = SSuper.block(0,     0,NTrial,NTrial);
        XTSigmaL = SSuper.block(NTrial,0,NTrial,NTrial);
      }
    } else {

    } 
  } // redDiag
  /** Form Residual Based Guess **/
  template<>
  void QuasiNewton<double>::formResidualGuess(double Omega, 
                                     const RealCMMap & ResR, RealCMMap & QR, 
                                     const RealCMMap & ResL, RealCMMap & QL){
    if(this->symmetrizedTrial_) {
      for(auto i = 0; i < this->N_; i++){
        QR(i) = ResR(i) * (*this->diag_)(i,0);
        QL(i) = ResL(i) * (*this->diag_)(i,0);
      }
      for(auto i = 0; i < this->N_/2; i++){
        QR(i) += Omega*ResL(i);
        QR(this->N_/2 + i) -= Omega*ResL(this->N_/2 + i);
        QL(i) += Omega*ResR(i);
        QL(this->N_/2 + i) -= Omega*ResR(this->N_/2 + i);
      }
      for(auto i = 0; i < this->N_; i++){
        QR(i) = QR(i) / (std::pow((*this->diag_)(i,0),2.0)-std::pow(Omega,2.0));
        QL(i) = QL(i) / (std::pow((*this->diag_)(i,0),2.0)-std::pow(Omega,2.0));
      }
    } else {
      for(auto i = 0; i < this->N_; i++)
        QR(i) = -ResR(i) / ((*this->diag_)(i,0) - Omega);
    }
  }
    /** Form New Perturbed Guess **/
  template<>
  void QuasiNewton<double>::formNewGuess(std::vector<bool> &resConv,int &NTrial, 
                                    int NNotConv, int &NOld, int &NNew){

    RealCMMap TrialVecR(this->TVecRMem,0,0);
    RealCMMap TrialVecL(this->TVecLMem,0,0);
    RealCMMap ResR     (this->ResRMem, 0,0);
    RealCMMap ResL     (this->ResLMem, 0,0);
    RealCMMap QR       (this->TVecRMem,0,0);
    RealCMMap RR       (this->ResRMem ,0,0);
    RealCMMap QL       (this->TVecLMem,0,0);
    RealCMMap RL       (this->ResLMem ,0,0);

    new (&TrialVecR) RealCMMap(this->TVecRMem,this->N_,NTrial+NNotConv);
    new (&ResR) RealCMMap(this->ResRMem,this->N_,NTrial);
    if(this->symmetrizedTrial_ || !this->isHermetian_){
      new (&TrialVecL) RealCMMap(this->TVecLMem,this->N_,NTrial+NNotConv);
      new (&ResL) RealCMMap(this->ResLMem,this->N_,NTrial);
    }

    RealVecMap ER(this->ERMem,NTrial);
    int INDX = 0;
    for(auto k = 0; k < this->nSek_; k++) {
      // If the residual for root "k" is not converged, construct
      // a perturbed guess vector according to Davidson's diagonal
      // preconditioning scheme.
      //
      // **WARNING** This is only valid for strongly diagonal dominant
      //             matricies. Convergence will be slow (maybe infinitely)
      //             if this criteria is not met.
      if(!resConv[k]) {
        new (&RR) RealCMMap(this->ResRMem + k*this->N_,this->N_,1);
        new (&QR) RealCMMap(this->TVecRMem+(NTrial+INDX)*this->N_,this->N_,1);
        if(this->symmetrizedTrial_ || !this->isHermetian_){
          new (&RL) RealCMMap(this->ResLMem + k*this->N_,this->N_,1);
          new (&QL) RealCMMap(this->TVecLMem+(NTrial+INDX)*this->N_,this->N_,1);
        }
         
        this->formResidualGuess(ER(k),RR,QR,RL,QL);
        INDX++;
      }
    }
    // Normalize and orthogonalize the new guess vectors to the
    // existing set using QR factorization
    int N,M,LDA,INFO;
    double * AMATR, * AMATL;
    N = TrialVecR.cols();
    M = TrialVecR.rows();
    LDA = TrialVecR.rows();
    AMATR = TrialVecR.data();
    if(this->symmetrizedTrial_ || !this->isHermetian_) AMATL = TrialVecL.data();
    double *TAU = this->LAPACK_SCR;
    this->WORK = TAU + N;
  
    dgeqrf_(&M,&N,AMATR,&LDA,TAU,this->WORK,&this->LWORK,&INFO);
    dorgqr_(&M,&N,&N,AMATR,&LDA,TAU,this->WORK,&this->LWORK,&INFO);
    if(this->symmetrizedTrial_ || !this->isHermetian_){
      dgeqrf_(&M,&N,AMATL,&LDA,TAU,this->WORK,&this->LWORK,&INFO);
      dorgqr_(&M,&N,&N,AMATL,&LDA,TAU,this->WORK,&this->LWORK,&INFO);
    }
    // Update number of vectors
    NOld = NTrial;
    NNew = NNotConv;
    NTrial += NNew;
  }
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

    // Symmetrize the trial vectors viz Kauczor et al. JCTC 7 (2010)
    TrialVecR = (*this->guessR_);
    if(this->symmetrizedTrial_){
      if(this->doRestart_) TrialVecL = (*this->guessL_);
      else                 TrialVecL = (*this->guessR_);
      if(!this->doRestart_){
        TrialVecR.block(this->N_/2,0,this->N_/2,this->nGuess_)
          = TrialVecR.block(0,0,this->N_/2,this->nGuess_);
        TrialVecL.block(this->N_/2,0,this->N_/2,this->nGuess_)
          = -TrialVecL.block(0,0,this->N_/2,this->nGuess_);
        // Normalize
        TrialVecR *= std::sqrt(0.5);
        TrialVecL *= std::sqrt(0.5);
      }
    }
    this->guessR_.reset();
    if(this->doRestart_) this->guessL_.reset();

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
//    if(this->doRestart_) CErr("QuasiNewton tried to extend the subspace");
      if(this->doRestart_) {
        this->allocGuess();
        RealCMMap UR(this->URMem,this->N_,this->nGuess_);
        (*this->guessR_) = UR;
        if(this->symmetrizedTrial_){
          RealCMMap UL(this->ULMem,this->N_,this->nGuess_);
          (*this->guessL_) = UL;
        } 
        std::memset(this->SCR,0.0,this->LenScr*sizeof(double));
        int N = this->guessR_->cols();
        int M = this->guessR_->rows();
        int LDA = this->guessR_->rows();
        int INFO;
        double *AMATR = this->guessR_->data();
        double *AMATL;
        if(this->symmetrizedTrial_ || !this->isHermetian_) AMATL = this->guessL_->data();
        double *TAU = this->LAPACK_SCR;
        this->WORK = TAU + N;
    
        dgeqrf_(&M,&N,AMATR,&LDA,TAU,this->WORK,&this->LWORK,&INFO);
        dorgqr_(&M,&N,&N,AMATR,&LDA,TAU,this->WORK,&this->LWORK,&INFO);
        if(this->symmetrizedTrial_ || !this->isHermetian_){
          dgeqrf_(&M,&N,AMATL,&LDA,TAU,this->WORK,&this->LWORK,&INFO);
          dorgqr_(&M,&N,&N,AMATL,&LDA,TAU,this->WORK,&this->LWORK,&INFO);
        }
//      CErr();
//      this->doRestart_ = false;
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
