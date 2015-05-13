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
void Davidson<RealMatrix>::runMicro(ostream &output ) {
  RealMatrix AX(this->n_,this->nGuess_);
  RealMatrix XTAX(this->nGuess_,this->nGuess_);
  RealMatrix U(this->n_,this->nGuess_);
  RealMatrix Res(this->n_,this->nGuess_);
  RealMatrix TrialVec(this->n_,this->nGuess_);
  RealMatrix T(this->n_,1);
  double *LAPACK_SCR;
  int LEN_LAPACK_SCR,LWORK;
  Eigen::Map<Eigen::VectorXd> ER(LAPACK_SCR,0);

  Eigen::SelfAdjointEigenSolver<RealMatrix> subDiag_;
  if(this->useLAPACK_) {
    LWORK = 6*this->n_;
    LEN_LAPACK_SCR = 0;
    LEN_LAPACK_SCR += this->maxSubSpace_; // Subspace eigenvalues (real)
//  LEN_LAPACK_SCR += this->maxSubSpace_*this->maxSubSpace_; // Subspace eigenvectors (right) (not needed for SY)
    LEN_LAPACK_SCR += LWORK; // LAPACK workspace
    LAPACK_SCR = new double[LEN_LAPACK_SCR];
  }

  int NTrial = this->nGuess_;
  TrialVec = (*this->guess_);
  for(auto iter = 0; iter < this->maxIter_; iter++){
    std::chrono::high_resolution_clock::time_point start,finish;
    std::chrono::duration<double> elapsed;
    output << "Starting Davidson Micro Iteration " << iter + 1 << endl;
    start = std::chrono::high_resolution_clock::now();

    // Matrix Product (AX). Keep around for reuse in computing
    // the residual vector
    AX = (*this->mat_) * TrialVec;  
   
    // Full projection of A onto subspace
    XTAX = TrialVec.transpose()*AX; 

    // Diagonalize the subspace
    if(!this->useLAPACK_) subDiag_.compute(XTAX);
    else {
      new (&ER) Eigen::Map<Eigen::VectorXd>(LAPACK_SCR,NTrial);
//    RealMap VR(LAPACK_SCR+NTrial,NTrial*NTrial); // DSYEV will overwrite XTAX
      char JOBZ = 'V';
      char UPLO = 'L';
      int INFO;
      dsyev_(&JOBZ,&UPLO,&NTrial,XTAX.data(),&NTrial,ER.data(),
             LAPACK_SCR+NTrial,&LWORK,&INFO); 
      if(INFO!=0) CErr("DSYEV failed to converge in Davison Iterations",output);
      XTAX.transposeInPlace(); // Convert to RowMajor...
    }
   
    
    // Reconstruct approximate eigenvectors
    if(!this->useLAPACK_) U = TrialVec * subDiag_.eigenvectors();
    else U = TrialVec * XTAX;

    // Stash away current approximation of eigenvalues and eigenvectors (NSek)
    if(!this->useLAPACK_) {
      (*this->eigenvalues_) = subDiag_.eigenvalues().block(0,0,this->nSek_,1);
    } else {
      (*this->eigenvalues_) = ER.block(0,0,this->nSek_,1);
    }
    (*this->eigenvector_) = U.block(0,0,this->n_,this->nSek_);
    
    // Construct the residual vector 
    // R = A*U - U*E = (AX)*c - U*E
    if(!this->useLAPACK_) {
      Res = AX*subDiag_.eigenvectors() - U*subDiag_.eigenvalues().asDiagonal();
    } else Res = AX*XTAX - U*ER.asDiagonal();

    // Vector to store convergence info
    std::vector<bool> resConv;
    int NNotConv = 0;

    // Loop over NSek residual vectors. Decide from which residuals
    // will be made perturbed guess vectors
    output << "  Checking Residual Norms:" << endl;
    for(auto k = 0; k < this->nSek_; k++) {
      if(Res.col(k).norm() < 5e-6) resConv.push_back(true);
      else {
        resConv.push_back(false); NNotConv++;
      }
      if(resConv[k]) {
        output << "    Norm of Residual " << k+1 << " = " << std::scientific 
               << Res.col(k).norm() << " \t \t Root has converged" <<endl;
      } else {
        output << "    Norm of Residual " << k+1 << " = " << std::scientific 
               << Res.col(k).norm() << " \t \t Root has not converged" <<endl;
      }
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
    TrialVec.conservativeResize(this->n_,NTrial+NNotConv);
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
          if(!this->useLAPACK_) {
            T(i,0) = - Res.col(k)(i) / ((*this->mat_)(i,i) - subDiag_.eigenvalues()(k));
          } else {
            T(i,0) = - Res.col(k)(i) / ((*this->mat_)(i,i) - ER(k));
          }
        }
//      output << TrialVec.rows() <<" " <<TrialVec.cols() << endl;
//      output << T.rows() <<" " <<T.cols() << endl;
        TrialVec.block(0,NTrial + INDX,this->n_,1) = T;
        INDX++;
      }
    }
    // Normalize and orthogonalize the new guess vectors to the
    // existing set using QR factorization (piviting a problem??)
    Eigen::FullPivHouseholderQR<RealMatrix> QR(TrialVec);
    TrialVec = QR.matrixQ().block(0,0,this->n_,NTrial+NNotConv);
    TrialVec = TrialVec*QR.colsPermutation().transpose(); // permute the vectors back
    
    NTrial += NNotConv;
    finish = std::chrono::high_resolution_clock::now();
    elapsed = finish - start;
    output << "Davidson Micro Iteration took " << std::fixed 
           << elapsed.count() << " secs" << endl << endl;
//  output << QR.colsPermutation().toDenseMatrix() << endl;
//  output << endl << TrialVec << endl;
  } 
  delete [] LAPACK_SCR; // Cleanup scratch space
}
} // namespace ChronusQ
