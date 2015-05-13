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
  Eigen::SelfAdjointEigenSolver<RealMatrix> subDiag_;

  int NTrial = this->nGuess_;
  TrialVec = (*this->guess_);
  for(auto iter = 0; iter < this->maxIter_; iter++){
   
    // Matrix Product (AX). Keep around for reuse in computing
    // the residual vector
    AX = (*this->mat_) * TrialVec;  
   
    // Full projection of A onto subspace
    XTAX = TrialVec.transpose()*AX; 

    // Diagonalize the subspace
    subDiag_.compute(XTAX);         // Diagonalize subspace
    
    // Reconstruct approximate eigenvectors
    U = TrialVec * subDiag_.eigenvectors();

    // Stash away current approximation of eigenvalues and eigenvectors (NSek)
    (*this->eigenvalues_) = subDiag_.eigenvalues().block(0,0,this->nSek_,1);
    (*this->eigenvector_) = U.block(0,0,this->n_,this->nSek_);
    
    // Construct the residual vector 
    // R = A*U - U*E = (AX)*c - U*E
    Res = AX*subDiag_.eigenvectors() - U*subDiag_.eigenvalues().asDiagonal();

    // Vector to store convergence info
    std::vector<bool> resConv;
    int NNotConv = 0;

    // Loop over NSek residual vectors. Decide from which residuals
    // will be made perturbed guess vectors
    for(auto k = 0; k < this->nSek_; k++) {
      if(Res.col(k).norm() < 1e-6) resConv.push_back(true);
      else {
        resConv.push_back(false); NNotConv++;
      }
      output << "Norm of Residual " << k+1 << " = " << Res.col(k).norm() 
           << " " << resConv[k] <<endl;
    }

    output << *this->eigenvalues_ << endl << endl;

    // If NSek (lowest) residuals are converged, exit, else, create
    // perturbed guess vectors for next iteration
    this->converged_ = (NNotConv == 0);
    if(this->converged_) break;


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
          T(i,0) = - Res.col(k)(i) / ((*this->mat_)(i,i) - subDiag_.eigenvalues()(k));
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
    
    NTrial += NNotConv;
  } 
}
} // namespace ChronusQ
