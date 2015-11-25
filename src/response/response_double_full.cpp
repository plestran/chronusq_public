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
#include <response.h>

namespace ChronusQ{
template<>
void Response<double>::fullFOPPA(){
  // Form Needed MO Integrals (Requires INCORE Ints)
  this->mointegrals_->formIAJB(false);
  this->mointegrals_->formIJAB(false);

  for(auto iMat = 0; iMat != iMatIter_.size(); iMat++){

    this->nSingleDim_ = this->nMatDim_[iMat];
    // LAPACK Scratch for Full Diagonalization
    char JOBZ  = 'V';
    char JOBVL = 'N';
    char UPLO  = 'U';
    int INFO;
    int LWORK = 6*this->nSingleDim_;
    double * WORK  = new double[LWORK]; 
    double * COPY, * WI;
    if(!this->doTDA_) {
      COPY = new double[this->nSingleDim_*this->nSingleDim_]; 
      WI   = new double[this->nSingleDim_];
    }

    if(this->iMatIter_[iMat] == FULL) {
      if(this->Ref_ != SingleSlater<double>::TCS) {
        /*
         *  Build the A Matrix (and B if RPA)
         *
         *  A(ia,s;jb,s') = delta_{ia;jb}*delta_{s;s'}*(eps(a,s) - eps(i,s))
         *                  + (ia,s | jb,s') - delta_{s;s'}*(ij,s | ab,s)
         *
         */ 

        int nSingOff = this->nSingleDim_/2; 
        int iMetricScale = -1;
        if(this->iMeth_ == STAB) iMetricScale = 1;
        // Loop over and populate All Alpha Block of A Matrix
        // ** If RPA, also populate all Alpha Block of B Matrix **
        for(auto i = 0, ia = 0; i < this->nOA_; i++)
        for(auto a = 0, A = this->nOA_; a < this->nVA_; a++, A++, ia++)
        for(auto j = 0, jb = 0; j < this->nOA_; j++)
        for(auto b = 0, B = this->nOA_; b < this->nVA_; b++, B++, jb++) {
          this->transDen_[iMat](ia,jb) = this->mointegrals_->IAJB(i,a,j,b) -  
                                         this->mointegrals_->IJAB(i,j,a,b);
          if(ia == jb)
            this->transDen_[iMat](ia,jb) += (*this->singleSlater_->epsA())(A) - 
                                            (*this->singleSlater_->epsA())(i);
          if(!this->doTDA_) {
            // Copy A into the lower right corner
            this->transDen_[iMat](ia+nSingOff,jb+nSingOff) =
              iMetricScale * this->transDen_[iMat](ia,jb);

            // Top Right B Matrix
            this->transDen_[iMat](ia,jb+nSingOff) =
              this->mointegrals_->IAJB(i,a,j,b) - 
              this->mointegrals_->IAJB(i,b,j,a);

            // Bottom Left B Matrix
            this->transDen_[iMat](ia+nSingOff,jb) =
              iMetricScale * this->transDen_[iMat](ia,jb+nSingOff);
          } // not TDA
        } // All Alpha Loop

        // Loop over and populate Mixed Alpha-Beta Block of A Matrix
        // ** If RPA, also populate Mixed Alpha-Beta Block of B Matrix **
        for(auto i = 0, ia = this->nOAVA_; i < this->nOB_; i++)
        for(auto a = 0, A = this->nOB_; a < this->nVB_; a++, A++, ia++)
        for(auto j = 0, jb = 0; j < this->nOA_; j++)
        for(auto b = 0, B = this->nOA_; b < this->nVA_; b++, B++, jb++) {
          this->transDen_[iMat](ia,jb) = this->mointegrals_->IAJB(i,a,j,b); 
          this->transDen_[iMat](jb,ia) = this->mointegrals_->IAJB(i,a,j,b); 
          if(!this->doTDA_) {
            // Copy A into the lower right corner
            this->transDen_[iMat](nSingOff+ia,nSingOff+jb) = 
              iMetricScale * this->transDen_[iMat](ia,jb); 
            this->transDen_[iMat](nSingOff+jb,nSingOff+ia) = 
              iMetricScale * this->transDen_[iMat](jb,ia); 

            // Top Right B Matrix
            this->transDen_[iMat](ia,nSingOff+jb) = 
              this->transDen_[iMat](ia,jb); 
            this->transDen_[iMat](jb,nSingOff+ia) = 
              this->transDen_[iMat](jb,ia); 

            // Bottom Left B Matrix
            this->transDen_[iMat](nSingOff+ia,jb) = 
              iMetricScale * this->transDen_[iMat](ia,jb); 
            this->transDen_[iMat](nSingOff+jb,ia) = 
              iMetricScale * this->transDen_[iMat](jb,ia); 
          } // Not TDA
        } // Mixed Alpha-Beta Loop

        // Loop over and populate All Beta Block of A Matrix
        // ** If RPA, also populate all Beta Block of B Matrix **
        for(auto i = 0, ia = this->nOAVA_; i < this->nOB_; i++)
        for(auto a = 0, A = this->nOB_; a < this->nVB_; a++, A++, ia++)
        for(auto j = 0, jb = this->nOAVA_; j < this->nOB_; j++)
        for(auto b = 0, B = this->nOB_; b < this->nVB_; b++, B++, jb++) {
          this->transDen_[iMat](ia,jb) = this->mointegrals_->IAJB(i,a,j,b) -  
                                         this->mointegrals_->IJAB(i,j,a,b);
          if(ia == jb)
            this->transDen_[iMat](ia,jb) += (*this->singleSlater_->epsA())(A) - 
                                            (*this->singleSlater_->epsA())(i);
          if(!this->doTDA_) {
            // Copy A into the lower right corner
            this->transDen_[iMat](ia+nSingOff,jb+nSingOff) =
              iMetricScale * this->transDen_[iMat](ia,jb);

            // Top Right B Matrix
            this->transDen_[iMat](ia,jb+nSingOff) =
              this->mointegrals_->IAJB(i,a,j,b) - 
              this->mointegrals_->IAJB(i,b,j,a);

            // Bottom Left B Matrix
            this->transDen_[iMat](ia+nSingOff,jb) =
            iMetricScale * this->transDen_[iMat](ia,jb+nSingOff);
          } // not TDA
        } // All Beta Loop
      } // TCS Check
    } // this->iMatIter_[iMat] == FULL
    else if(this->iMatIter_[iMat] == SINGLETS){
      for(auto i = 0, ia = 0; i < this->nOA_; i++)
      for(auto a = 0, A = this->nOA_; a < this->nVA_; a++, A++, ia++)
      for(auto j = 0, jb = 0; j < this->nOA_; j++)
      for(auto b = 0, B = this->nOA_; b < this->nVA_; b++, B++, jb++) {
        this->transDen_[iMat](ia,jb) = 2*this->mointegrals_->IAJB(i,a,j,b) -  
                                         this->mointegrals_->IJAB(i,j,a,b);
        if(ia == jb)
          this->transDen_[iMat](ia,jb) += (*this->singleSlater_->epsA())(A) - 
                                          (*this->singleSlater_->epsA())(i);
      }
    } // Singlets
    else if(this->iMatIter_[iMat] == TRIPLETS){
      for(auto i = 0, ia = 0; i < this->nOA_; i++)
      for(auto a = 0, A = this->nOA_; a < this->nVA_; a++, A++, ia++)
      for(auto j = 0, jb = 0; j < this->nOA_; j++)
      for(auto b = 0, B = this->nOA_; b < this->nVA_; b++, B++, jb++) {
        this->transDen_[iMat](ia,jb) = - this->mointegrals_->IJAB(i,j,a,b);
        if(ia == jb)
          this->transDen_[iMat](ia,jb) += (*this->singleSlater_->epsA())(A) - 
                                          (*this->singleSlater_->epsA())(i);
      }
    } // Triplets

    int N = this->nSingleDim_;
    if(this->doTDA_ || this->iMeth_ == STAB)
      // Diagonalize A or ABBA (for stability)
      dsyev_(&JOBZ,&UPLO,&N,this->transDen_[iMat].data(),&N,
        this->frequencies_[iMat].data(),WORK,&LWORK,&INFO);
    else {
      // Copy the response Matrix to a Temporary for DGEEV
      std::memcpy(COPY,this->transDen_[iMat].data(),
        this->nSingleDim_*this->nSingleDim_*sizeof(double));

      // Zero out Eigenvalue Storage
      std::memset(WI,0.0,this->nSingleDim_*sizeof(double));

      // Solve full eigenvalue problem for FOPPA
      dgeev_(&JOBVL,&JOBZ,&N,COPY,&N,this->frequencies_[iMat].data(),WI,
        this->transDen_[iMat].data(),&N,this->transDen_[iMat].data(),&N,
        WORK,&LWORK,&INFO);

      // Sort the eigenvalues
      std::sort(this->frequencies_[iMat].data(),
                this->frequencies_[iMat].data()+
                this->frequencies_[iMat].size());

      // Cleanup DGEEV required memory
      delete [] COPY;
      delete [] WI;
    }

    // Check if the diagonalization Converged
    if(INFO != 0){
      std::string msg;
      msg = "Full In-Core Diagonalization ";
      if(this->iMeth_ != RPA)
        msg += "(DSYEV) ";
      else
        msg += "(DGEEV) ";
      msg += "Failed for IMat = "; 
      msg += std::to_string(static_cast<int>(this->iMatIter_[iMat]));
      CErr(msg,this->fileio_->out);
    };

    RealVecMap Eig(this->frequencies_[iMat].data(),this->nSingleDim_);
    prettyPrint(this->fileio_->out,Eig*phys.eVPerHartree,"Eig");
    // Cleanup LAPACK Memory
    delete[] WORK;

  } // loop over iMat 
}; // fullFOPPA (T = double)

template<>
void Response<double>::fullPPRPA(){
  // Form Needed MO Integrals (Requires INCORE Ints)
  this->mointegrals_->formABCD(false);
  this->mointegrals_->formIAJB(false);
  this->mointegrals_->formIJKL(false);

  for(auto iMat = 0; iMat != iMatIter_.size(); iMat++){

    this->nSingleDim_ = this->nMatDim_[iMat];
    // LAPACK Scratch for Full Diagonalization
    char JOBZ  = 'V';
    char JOBVL = 'N';
    char UPLO  = 'U';
    int INFO;
    int LWORK = 6*this->nSingleDim_;
    double * WORK  = new double[LWORK]; 
    double * COPY, * WI;
    if(!this->doTDA_) {
      COPY = new double[this->nSingleDim_*this->nSingleDim_]; 
      WI   = new double[this->nSingleDim_];
    }

    if(this->iMatIter_[iMat] == AA_PPRPA){ 
      auto iOff = this->nVAVA_SLT_;
      auto iMetricScale = -1;
      // Place A Matrix
      for(auto a = 0, ab = 0; a < this->nVA_; a++      )
      for(auto b = 0        ; b < a        ;  b++, ab++)
      for(auto c = 0, cd = 0; c < this->nVA_; c++      )
      for(auto d = 0        ; d < c        ;  d++, cd++){
        this->transDen_[iMat](ab,cd) = this->mointegrals_->ABCD(a,c,b,d) -
                                       this->mointegrals_->ABCD(a,d,b,c);  
        if(ab == cd) this->transDen_[iMat](ab,cd) +=
          (*this->singleSlater_->epsA())(a+this->nOA_) + 
          (*this->singleSlater_->epsA())(b+this->nOA_) - 2*this->rMu_; 
      }

      // Place C Matrix
      for(auto i = 0, ij = iOff; i < this->nOA_; i++      )
      for(auto j = 0           ; j < i         ; j++, ij++)
      for(auto k = 0, kl = iOff; k < this->nOA_; k++      )
      for(auto l = 0           ; l < k         ; l++, kl++){
        this->transDen_[iMat](ij,kl) = this->mointegrals_->IJKL(i,k,j,l) - 
                                       this->mointegrals_->IJKL(i,l,j,k);
        if(ij == kl) this->transDen_[iMat](ij,kl) -=
          (*this->singleSlater_->epsA())(i) + 
          (*this->singleSlater_->epsA())(j) - 2*this->rMu_; 

        this->transDen_[iMat](ij,kl) *= iMetricScale;
      }

      // Place B Matricies
      for(auto a = 0, ab =    0; a < this->nVA_; a++      )
      for(auto b = 0           ; b < a        ;  b++, ab++)
      for(auto i = 0, ij = iOff; i < this->nOA_; i++      )
      for(auto j = 0           ; j < i         ; j++, ij++) {
        this->transDen_[iMat](ab,ij) = this->mointegrals_->IAJB(i,a,j,b) - 
                                       this->mointegrals_->IAJB(i,b,j,a);
        this->transDen_[iMat](ij,ab) = 
          iMetricScale*this->transDen_[iMat](ab,ij);
      }
    } // this->iMatIter_[iMat] == AA_PPRPA
    else if(this->iMatIter_[iMat] == AB_PPRPA){ 
      auto iOff = this->nVAVB_;
      auto iMetricScale = -1;
     
      // Place A Matrix
      for(auto a = 0, ab = 0; a < this->nVA_; a++      )
      for(auto b = 0        ; b < this->nVB_; b++, ab++)
      for(auto c = 0, cd = 0; c < this->nVA_; c++      )
      for(auto d = 0        ; d < this->nVB_; d++, cd++){
        this->transDen_[iMat](ab,cd) = this->mointegrals_->ABCD(a,c,b,d);
        if(a == c && b == d) this->transDen_[iMat](ab,cd) +=
          (*this->singleSlater_->epsA())(a+this->nOA_) + 
          (*this->singleSlater_->epsA())(b+this->nOB_) - 2*this->rMu_; 
      }

      // Place C Matrix
      for(auto i = 0, ij = iOff; i < this->nOA_; i++      )
      for(auto j = 0           ; j < this->nOB_; j++, ij++)
      for(auto k = 0, kl = iOff; k < this->nOA_; k++      )
      for(auto l = 0           ; l < this->nOB_; l++, kl++){
        this->transDen_[iMat](ij,kl) = this->mointegrals_->IJKL(i,k,j,l); 

        if(i == k && j == l) this->transDen_[iMat](ij,kl) -=
          (*this->singleSlater_->epsA())(i) + 
          (*this->singleSlater_->epsA())(j) - 2*this->rMu_; 

        this->transDen_[iMat](ij,kl) *= iMetricScale;
      }

      // Place B Matricies
      for(auto a = 0, ab = 0; a < this->nVA_; a++      )
      for(auto b = 0        ; b < this->nVB_; b++, ab++)
      for(auto i = 0, ij = iOff; i < this->nOA_; i++      )
      for(auto j = 0           ; j < this->nOB_; j++, ij++){
        this->transDen_[iMat](ab,ij) = this->mointegrals_->IAJB(i,a,j,b); 
        this->transDen_[iMat](ij,ab) = 
          iMetricScale*this->transDen_[iMat](ab,ij);
      }
    } // this->iMatIter_[iMat] == AB_PPRPA
    else if(this->iMatIter_[iMat] == BB_PPRPA){ 
    } // this->iMatIter_[iMat] == BB_PPRPA
    else if(this->iMatIter_[iMat] == AAA_PPTDA){
      for(auto a = 0, ab = 0; a < this->nVA_; a++      )
      for(auto b = 0        ; b < a        ;  b++, ab++)
      for(auto c = 0, cd = 0; c < this->nVA_; c++      )
      for(auto d = 0        ; d < c        ;  d++, cd++){
        this->transDen_[iMat](ab,cd) = this->mointegrals_->ABCD(a,c,b,d) -
                                       this->mointegrals_->ABCD(a,d,b,c);  
        if(ab == cd) this->transDen_[iMat](ab,cd) +=
          (*this->singleSlater_->epsA())(a+this->nOA_) + 
          (*this->singleSlater_->epsA())(b+this->nOA_) - 2*this->rMu_; 
      }
    } // this->iMatIter_[iMat] == AAA_PPTDA
    else if(this->iMatIter_[iMat] == AAB_PPTDA){
      for(auto a = 0, ab = 0; a < this->nVA_; a++      )
      for(auto b = 0        ; b < this->nVB_; b++, ab++)
      for(auto c = 0, cd = 0; c < this->nVA_; c++      )
      for(auto d = 0        ; d < this->nVB_; d++, cd++){
        this->transDen_[iMat](ab,cd) = this->mointegrals_->ABCD(a,c,b,d);
        if(a == c && b == d) this->transDen_[iMat](ab,cd) +=
          (*this->singleSlater_->epsA())(a+this->nOA_) + 
          (*this->singleSlater_->epsA())(b+this->nOB_) - 2*this->rMu_; 
      }
    } // this->iMatIter_[iMat] == AAB_PPTDA
    else if(this->iMatIter_[iMat] == ABB_PPTDA){
    } // this->iMatIter_[iMat] == ABB_PPTDA
    else if(this->iMatIter_[iMat] == CAA_PPTDA){
      for(auto i = 0, ij = 0; i < this->nOA_; i++      )
      for(auto j = 0        ; j < i         ; j++, ij++)
      for(auto k = 0, kl = 0; k < this->nOA_; k++      )
      for(auto l = 0        ; l < k         ; l++, kl++){
        this->transDen_[iMat](ij,kl) = this->mointegrals_->IJKL(i,k,j,l) - 
                                       this->mointegrals_->IJKL(i,l,j,k);
        if(ij == kl) this->transDen_[iMat](ij,kl) -=
          (*this->singleSlater_->epsA())(i) + 
          (*this->singleSlater_->epsA())(j) - 2*this->rMu_; 
      }
    } // this->iMatIter_[iMat] == CAA_PPTDA
    else if(this->iMatIter_[iMat] == CAB_PPTDA){
      for(auto i = 0, ij = 0; i < this->nOA_; i++      )
      for(auto j = 0        ; j < this->nOB_; j++, ij++)
      for(auto k = 0, kl = 0; k < this->nOA_; k++      )
      for(auto l = 0        ; l < this->nOB_; l++, kl++){
        this->transDen_[iMat](ij,kl) = this->mointegrals_->IJKL(i,k,j,l); 
        if(i == k && j == l) this->transDen_[iMat](ij,kl) -=
          (*this->singleSlater_->epsA())(i) + 
          (*this->singleSlater_->epsA())(j) - 2*this->rMu_; 
      }
    } // this->iMatIter_[iMat] == CAB_PPTDA
    else if(this->iMatIter_[iMat] == CBB_PPTDA){
    } // this->iMatIter_[iMat] == CBB_PPTDA

    int N = this->nSingleDim_;
    if(this->doTDA_)
      dsyev_(&JOBZ,&UPLO,&N,this->transDen_[iMat].data(),&N,
        this->frequencies_[iMat].data(),WORK,&LWORK,&INFO);
    else {
      // Copy the response Matrix to a Temporary for DGEEV
      std::memcpy(COPY,this->transDen_[iMat].data(),
        this->nSingleDim_*this->nSingleDim_*sizeof(double));

      // Zero out Eigenvalue Storage
      std::memset(WI,0.0,this->nSingleDim_*sizeof(double));

      // Solve full eigenvalue problem for FOPPA
      dgeev_(&JOBVL,&JOBZ,&N,COPY,&N,this->frequencies_[iMat].data(),WI,
        this->transDen_[iMat].data(),&N,this->transDen_[iMat].data(),&N,
        WORK,&LWORK,&INFO);

      this->frequencies_[iMat] *= -1;
      // Sort the eigenvalues
      std::sort(this->frequencies_[iMat].data(),
                this->frequencies_[iMat].data()+
                this->frequencies_[iMat].size());
      this->frequencies_[iMat] *= -1;

      // Cleanup DGEEV required memory
      delete [] COPY;
      delete [] WI;
    }

    // Check if the diagonalization Converged
    if(INFO != 0){
      std::string msg;
      msg = "Full In-Core Diagonalization ";
      if(this->iMeth_ != PPRPA)
        msg += "(DSYEV) ";
      else
        msg += "(DGEEV) ";
      msg += "Failed for IMat = "; 
      msg += std::to_string(static_cast<int>(this->iMatIter_[iMat]));
      CErr(msg,this->fileio_->out);
    };

    RealVecMap Eig(this->frequencies_[iMat].data(),this->nSingleDim_);
    prettyPrint(this->fileio_->out,Eig,"Eig");
    // Cleanup LAPACK Memory
    delete[] WORK;
  }; // loop over iMat
}; // fullPPRPA

}; // namespace ChronusQ
