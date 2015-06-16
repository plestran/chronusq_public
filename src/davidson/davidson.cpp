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
#include <davidson.h>

namespace ChronusQ {
template<>
void Davidson<double>::runMicro(ostream &output ) {

  /* Explicit Memory Allocation */

  int LenScr  = 0;
  int LenSigma   = this->n_ * this->maxSubSpace_;
  int LenXTSigma = this->maxSubSpace_ * this->maxSubSpace_;
  int LenU    = this->n_ * this->maxSubSpace_;
  int LenRes  = this->n_ * this->maxSubSpace_;
  int LenTVec = this->n_ * this->maxSubSpace_;
  int LenT    = this->n_;
  int LEN_LAPACK_SCR = 0;
  int LWORK = 0;

  /* *******************************************************
   * *********Determine length of scratch space*************
   * *******************************************************
   *
   * (1) Linear transform of A onto right / gerade trial vectors 
   *
   * (2) Reduced dimension of A onto the right / gerade trial vector subspace
   *     (also stores the reduced dimension eigenvectors after diagonalization
   *      via DSYEV)
   *
   * (3) Approximate right / gerade eigenvectors 
   *
   * (4) Residuals of the right / gerade eigenvectors 
   *
   * (5) Right / gerade trial vectors 
   *
   * (6) Temp storage for a single precondictioned eigenvector 
   *
   * (12) Local copy of the real part of the eigenvalues (reused for Tau storage for QR)
   *
   * (13) Length of LAPACK workspace (used in all LAPACK Calls)
   *
   * (14) Total double precision words required for LAPACK
   */

  LenScr += LenSigma;   // 1
  LenScr += LenXTSigma; // 2
  LenScr += LenU;       // 3
  LenScr += LenRes;     // 4 
  LenScr += LenTVec;    // 5
  LenScr += LenT;       // 6


  if(!this->hermetian_ || this->symmetrized_){
    LenScr += LenSigma;   // 7
    LenScr += LenXTSigma; // 8
    LenScr += LenU;       // 9
    LenScr += LenRes;     // 10
    LenScr += LenTVec;    // 11
  }

  LWORK = 6*this->n_;
  LEN_LAPACK_SCR += this->maxSubSpace_; // 12
  LEN_LAPACK_SCR += LWORK;              // 13
  LenScr += LEN_LAPACK_SCR;             // 14

  double * SCR         = NULL; 
  double * SigmaRMem   = NULL; 
  double * SigmaLMem   = NULL; 
  double * XTSigmaRMem = NULL; 
  double * XTSigmaLMem = NULL; 
  double * RhoRMem     = NULL; 
  double * RhoLMem     = NULL; 
  double * XTRhoRMem   = NULL; 
  double * XTRhoLMem   = NULL; 
  double * URMem       = NULL; 
  double * ULMem       = NULL; 
  double * ResRMem     = NULL; 
  double * ResLMem     = NULL; 
  double * TVecRMem    = NULL;         
  double * TVecLMem    = NULL;         
  double * TMem        = NULL;        
  double * LAPACK_SCR  = NULL;

  // Allocate Scratch Space
  SCR = new double [LenScr]; 

  // Partition scratch space
  SigmaRMem     = SCR;
  XTSigmaRMem   = SigmaRMem   + LenSigma;
  URMem         = XTSigmaRMem + LenXTSigma;
  ResRMem       = URMem       + LenU; 
  TVecRMem      = ResRMem     + LenRes;
  TMem          = TVecRMem    + LenTVec;
  LAPACK_SCR    = TMem        + LenT;

//RealMatrix TrialVecR(this->n_,this->nGuess_);

  RealCMMap SigmaR(   SigmaRMem,  0,0);
  RealCMMap NewSR(    SigmaRMem,  0,0);
  RealCMMap XTSigmaR( XTSigmaRMem,0,0);
  RealCMMap RhoR(     RhoRMem,    0,0);
  RealCMMap NewRhoR(  RhoRMem,    0,0);
  RealCMMap UR(       URMem,      0,0);
  RealCMMap ResR(     ResRMem,    0,0);
  RealCMMap TrialVecR(TVecRMem,   0,0);
  RealCMMap NewVecR(  TVecRMem,   0,0);

  RealCMMap SigmaL(   SigmaLMem,  0,0);
  RealCMMap NewSL(    SigmaLMem,  0,0);
  RealCMMap XTSigmaL( XTSigmaLMem,0,0);
  RealCMMap RhoL(     RhoLMem,    0,0);
  RealCMMap NewRhoL(  RhoLMem,    0,0);
  RealCMMap UL(       ULMem,      0,0);
  RealCMMap ResL(     ResLMem,    0,0);
  RealCMMap TrialVecL(TVecLMem,   0,0);
  RealCMMap NewVecL(  TVecLMem,   0,0);

  RealCMMap T(TMem,this->n_,  1);


  RealVecMap ER(LAPACK_SCR,0);
  RealVecMap EI(LAPACK_SCR,0);
  RealMap    VR(LAPACK_SCR,0,0);
  RealMap    VL(LAPACK_SCR,0,0);
/*
  if(!this->hermetian_) {
    // Subspace eigenvalues (imag)
    LEN_LAPACK_SCR += this->maxSubSpace_; 
    // Subspace eigenvectors (right) (not needed for SY)
    LEN_LAPACK_SCR += this->maxSubSpace_*this->maxSubSpace_; 
    // Subspace eigenvectors (left) (not needed for SY)
    LEN_LAPACK_SCR += this->maxSubSpace_*this->maxSubSpace_;
  }
*/
  char JOBVR = 'V';
  char JOBVL = 'N';
  char UPLO = 'L';
  int INFO;

  int NTrial = this->nGuess_;
  int NOld   = 0;
  int NNew   = this->nGuess_;
  new (&TrialVecR) RealCMMap(TVecRMem,this->n_,NTrial);

  TrialVecR = (*this->guess_);

  for(auto iter = 0; iter < this->maxIter_; iter++){
    std::chrono::high_resolution_clock::time_point start,finish;
    std::chrono::duration<double> elapsed;
    output << "Starting Davidson Micro Iteration " << iter + 1 << endl;
    start = std::chrono::high_resolution_clock::now();

    // Resize the Eigen Maps to fit new vectors
    new (&SigmaR)   RealCMMap(SigmaRMem,  this->n_,NTrial);
    new (&XTSigmaR) RealCMMap(XTSigmaRMem,NTrial,  NTrial);
    new (&UR)       RealCMMap(URMem,      this->n_,NTrial);
    new (&ResR)     RealCMMap(ResRMem,    this->n_,NTrial);
    new (&NewSR)    RealCMMap(SigmaRMem+NOld*this->n_,this->n_,NNew);
    new (&NewVecR)  RealCMMap(TVecRMem+NOld*this->n_,this->n_,NNew);

    // Matrix Product (Sigma / AX). Keep around for reuse in computing
    // the residual vector
    if(this->method_ == SDResponse::CIS) 
      this->sdr_->formRM3(NewVecR,NewSR,NewRhoR);
    else if(this->AX_==NULL) NewSR = (*this->mat_) * NewVecR;  
    else NewSR = this->AX_(*this->mat_,NewVecR);
//  cout << SigmaR << endl << endl;
   
    // Full projection of A onto subspace
    XTSigmaR = TrialVecR.transpose()*SigmaR; 

    // Diagonalize the subspace
    new (&ER) RealVecMap(LAPACK_SCR,NTrial);
    int IOff = NTrial;
/*
    if(!this->hermetian_) {
      new (&EI) Eigen::Map<Eigen::VectorXd>(LAPACK_SCR+IOff,NTrial);
      IOff += NTrial;
      new (&VR) RealMap(LAPACK_SCR+IOff,NTrial,NTrial);
      IOff += NTrial*NTrial;
    }
*/
    if(this->hermetian_) {
      dsyev_(&JOBVR,&UPLO,&NTrial,XTSigmaR.data(),&NTrial,ER.data(),
             LAPACK_SCR+IOff,&LWORK,&INFO); 
      if(INFO!=0) CErr("DSYEV failed to converge in Davison Iterations",output);
    } 
/*
    else {
      dgeev_(&JOBVL,&JOBVR,&NTrial,XTSigmaR.data(),&NTrial,ER.data(),
             EI.data(),VL.data(),&NTrial,VR.data(),&NTrial,
             LAPACK_SCR+IOff,&LWORK,&INFO);
      if(INFO!=0) CErr("DGEEV failed to converge in Davison Iterations",output);
      VR.transposeInPlace(); // Convert ro RowMajor...
//    eigSort(false,ER,VR,VL);
    }
*/
    
   
    
    // Reconstruct approximate eigenvectors
    if(this->hermetian_) UR = TrialVecR * XTSigmaR;
//  else                 UR = TrialVecR * VR;

    // Stash away current approximation of eigenvalues and eigenvectors (NSek)
    (*this->eigenvalues_) = ER.block(0,0,this->nSek_,1);
    (*this->eigenvector_) = UR.block(0,0,this->n_,this->nSek_);

    
    // Construct the residual vector 
    // R = A*U - U*E = (AX)*c - U*E
    if(this->hermetian_) ResR = SigmaR*XTSigmaR - UR*ER.asDiagonal();
/*
    else {
      int NImag = 0;
      for(auto k = 0; k < NTrial; k++) if(std::abs(EI(k))<1e-10) NImag++;
      ResR.resize(this->n_,NTrial+NImag);
      for(auto k = 0; k < NTrial+NImag; k++) {
        ResR.col(k) = SigmaR*VR.col(k) - ER(k)*UR.col(k);
        if(std::abs(EI(k))<1e-10) {
          ResR.col(k+1) = SigmaR*VR.col(k) - EI(k)*UR.col(k);
          k++;
        }
      }
    }
*/
    

    // Vector to store convergence info
    std::vector<bool> resConv;
    int NNotConv = 0;

    // Loop over NSek residual vectors. Decide from which residuals
    // will be made perturbed guess vectors
    for(auto k = 0; k < this->nSek_; k++) {
      if(ResR.col(k).norm() < 5e-6) resConv.push_back(true);
      else {
        resConv.push_back(false); NNotConv++;
      }
    }
    output << "  Checking Davidson Convergence:" << endl;
    output << "    " << std::setw(8)  << " " << std::setw(32) << std::left << "    Roots at Current Iteration:";
    output << std::setw(32) << std::left << "    Norm of Residual:" << endl;
    for(auto k = 0 ; k < this->nSek_; k++){
      output << "    " << std::setw(12) << "State " + std::to_string(k+1) + ":";
      output << std::setw(32) << std::left << std::fixed << (*this->eigenvalues_)(k,0);
      output << std::setw(15) << std::left << std::scientific << ResR.col(k).norm();
      if(resConv[k]) output << "     Root has converged" << endl;
      else output << "     Root has not converged" << endl;
    }

//  output << *this->eigenvalues_ << endl << endl;

    // If NSek (lowest) residuals are converged, exit, else, create
    // perturbed guess vectors for next iteration
    this->converged_ = (NNotConv == 0);
    if(this->converged_) {
      finish = std::chrono::high_resolution_clock::now();
      elapsed = finish - start;
      output << "  Davidson Micro Iteration took " << std::fixed 
             << elapsed.count() << " secs" << endl << endl;
      break;
    }

    // Resize the trial vector dimension to contain the new perturbed
    // guess vectors
    new (&TrialVecR) RealCMMap(TVecRMem,this->n_,NTrial+NNotConv);
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
        for(auto i = 0; i < this->n_; i++) {
          if(this->method_ == SDResponse::CIS) {
            T(i,0) = - ResR.col(k)(i) / ((*this->sdr_->rmDiag())(i,0) - ER(k));
          } else {
            T(i,0) = - ResR.col(k)(i) / ((*this->mat_)(i,i) - ER(k));
          }
        }
        TrialVecR.block(0,NTrial + INDX,this->n_,1) = T;
        INDX++;
      }
    }
    // Normalize and orthogonalize the new guess vectors to the
    // existing set using QR factorization
    int N = TrialVecR.cols();
    int M = TrialVecR.rows();
    dgeqrf_(&M,&N,TrialVecR.data(),&M,LAPACK_SCR,LAPACK_SCR+N,&LWORK,&INFO);
    dorgqr_(&M,&N,&N,TrialVecR.data(),&M,LAPACK_SCR,LAPACK_SCR+N,&LWORK,&INFO);
    
    NOld = NTrial;
    NNew = NNotConv;
    NTrial += NNotConv;
    finish = std::chrono::high_resolution_clock::now();
    elapsed = finish - start;
    output << "Davidson Micro Iteration took " << std::fixed 
           << elapsed.count() << " secs" << endl << endl;
  } 
  delete [] SCR; // Cleanup scratch space
}
} // namespace ChronusQ
