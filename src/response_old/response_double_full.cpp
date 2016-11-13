/* 
 *  This file is part of the Chronus Quantum (ChronusQ) software package
 *  
 *  Copyright (C) 2014-2016 Li Research Group (University of Washington)
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
    this->currentMat_ = this->iMatIter_[iMat];
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

    int nSingOff = this->nSingleDim_/2; 
    int iMetricScale = -1;
    if(this->currentMat_ == FULL) {
      if(this->Ref_ != SingleSlater<double>::TCS) {
        /*
         *  Build the A Matrix (and B if RPA)
         *
         *  A(ia,s;jb,s') = delta_{ia;jb}*delta_{s;s'}*(eps(a,s) - eps(i,s))
         *                  + (ia,s | jb,s') - delta_{s;s'}*(ij,s | ab,s)
         *
         */ 

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
    } // this->currentMat_ == FULL
    else if(this->currentMat_ == SINGLETS){
      for(auto i = 0, ia = 0; i < this->nOA_; i++)
      for(auto a = 0, A = this->nOA_; a < this->nVA_; a++, A++, ia++)
      for(auto j = 0, jb = 0; j < this->nOA_; j++)
      for(auto b = 0, B = this->nOA_; b < this->nVA_; b++, B++, jb++) {
        this->transDen_[iMat](ia,jb) = 2*this->mointegrals_->IAJB(i,a,j,b) -  
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
            2*this->mointegrals_->IAJB(i,a,j,b) - 
            this->mointegrals_->IAJB(i,b,j,a);
       
          // Bottom Left B Matrix
          this->transDen_[iMat](ia+nSingOff,jb) =
            iMetricScale * this->transDen_[iMat](ia,jb+nSingOff);
        } // not TDA
      }
    } // Singlets
    else if(this->currentMat_ == TRIPLETS){
      for(auto i = 0, ia = 0; i < this->nOA_; i++)
      for(auto a = 0, A = this->nOA_; a < this->nVA_; a++, A++, ia++)
      for(auto j = 0, jb = 0; j < this->nOA_; j++)
      for(auto b = 0, B = this->nOA_; b < this->nVA_; b++, B++, jb++) {
        this->transDen_[iMat](ia,jb) = - this->mointegrals_->IJAB(i,j,a,b);
        if(ia == jb)
          this->transDen_[iMat](ia,jb) += (*this->singleSlater_->epsA())(A) - 
                                          (*this->singleSlater_->epsA())(i);
        if(!this->doTDA_) {
          // Copy A into the lower right corner
          this->transDen_[iMat](ia+nSingOff,jb+nSingOff) =
            iMetricScale * this->transDen_[iMat](ia,jb);
       
          // Top Right B Matrix
          this->transDen_[iMat](ia,jb+nSingOff) =
            -this->mointegrals_->IAJB(i,b,j,a);
       
          // Bottom Left B Matrix
          this->transDen_[iMat](ia+nSingOff,jb) =
            iMetricScale * this->transDen_[iMat](ia,jb+nSingOff);
        } // not TDA
      }
    } // Triplets

    RealMatrix CPY(this->transDen_[iMat]);
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
      msg += std::to_string(static_cast<int>(this->currentMat_));
      CErr(msg,this->fileio_->out);
    };

  //RealVecMap Eig(this->frequencies_[iMat].data(),this->nSingleDim_);
  //prettyPrint(this->fileio_->out,Eig*phys.eVPerHartree,"Eig");
    // Cleanup LAPACK Memory
    delete[] WORK;

    this->nSek_ = this->nSingleDim_;
  //RealMap T(this->transDen_[iMat].data(),this->nSingleDim_,this->nSek_);
  //RealMatrix Sigma(this->nSingleDim_,2*this->nSek_);
  //RealMap SRVec(Sigma.data(),this->nSingleDim_,this->nSek_);
  //RealMap SLVec(Sigma.data()+this->nSek_*this->nSingleDim_,this->nSingleDim_,this->nSek_);
  //RealMatrix Rho(this->nSingleDim_,2*this->nSek_);
  //RealMap RRVec(Rho.data(),this->nSingleDim_,this->nSek_);
  //RealMap RLVec(Rho.data()+this->nSek_*this->nSingleDim_,this->nSingleDim_,this->nSek_);

  //prettyPrint(cout,T,"T");
  //this->linearTransFOPPA(T,T,SRVec,SLVec,RRVec,RLVec);

  //CPY.block(this->nSingleDim_/2,0,this->nSingleDim_/2,this->nSingleDim_/2) *= -1;
  //CPY.block(this->nSingleDim_/2,this->nSingleDim_/2,this->nSingleDim_/2,this->nSingleDim_/2) *= -1;
  //prettyPrint(cout,CPY*T,"Correct");
  //prettyPrint(cout,SRVec,"Test S1");
  //prettyPrint(cout,SLVec,"Test S2");
  //prettyPrint(cout,RRVec,"Test R1");
  //prettyPrint(cout,RLVec,"Test R2");
  //prettyPrint(cout,CPY*T-SRVec,"DIFF 1");
  //prettyPrint(cout,CPY*T-SLVec,"DIFF 2");
  //prettyPrint(cout,RRVec - RLVec,"DIFF 3");

  } // loop over iMat 
}; // fullFOPPA (T = double)

template<>
void Response<double>::fullPPRPA(){
  // Form Needed MO Integrals (Requires INCORE Ints)
  this->mointegrals_->formABCD(false);
  if(!this->doTDA_) this->mointegrals_->formIAJB(false);
  this->mointegrals_->formIJKL(false);

  this->rMu_ = 0.0;
  for(auto iMat = 0; iMat != iMatIter_.size(); iMat++){

    this->currentMat_ = this->iMatIter_[iMat];
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

    if(this->currentMat_ == AA_PPRPA){ 
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
    } // this->currentMat_ == AA_PPRPA
    else if(this->currentMat_ == AB_PPRPA){ 
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
    } // this->currentMat_ == AB_PPRPA
    else if(this->currentMat_ == BB_PPRPA){ 
    // ************************** //
    // ** BETA-BETA BLOCKS NYI ** //
    // ************************** //
    } // this->currentMat_ == BB_PPRPA
    else if(this->currentMat_ == AAA_PPTDA){
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
    } // this->currentMat_ == AAA_PPTDA
    else if(this->currentMat_ == AAB_PPTDA){
      for(auto a = 0, ab = 0; a < this->nVA_; a++      )
      for(auto b = 0        ; b < this->nVB_; b++, ab++)
      for(auto c = 0, cd = 0; c < this->nVA_; c++      )
      for(auto d = 0        ; d < this->nVB_; d++, cd++){
        this->transDen_[iMat](ab,cd) = this->mointegrals_->ABCD(a,c,b,d);
        if(a == c && b == d) this->transDen_[iMat](ab,cd) +=
          (*this->singleSlater_->epsA())(a+this->nOA_) + 
          (*this->singleSlater_->epsA())(b+this->nOB_) - 2*this->rMu_; 
      }
    } // this->currentMat_ == AAB_PPTDA
    else if(this->currentMat_ == ABB_PPTDA){
    // ************************** //
    // ** BETA-BETA BLOCKS NYI ** //
    // ************************** //
    } // this->currentMat_ == ABB_PPTDA
    else if(this->currentMat_ == CAA_PPTDA){
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
    } // this->currentMat_ == CAA_PPTDA
    else if(this->currentMat_ == CAB_PPTDA){
      for(auto i = 0, ij = 0; i < this->nOA_; i++      )
      for(auto j = 0        ; j < this->nOB_; j++, ij++)
      for(auto k = 0, kl = 0; k < this->nOA_; k++      )
      for(auto l = 0        ; l < this->nOB_; l++, kl++){
        this->transDen_[iMat](ij,kl) = this->mointegrals_->IJKL(i,k,j,l); 
        if(i == k && j == l) this->transDen_[iMat](ij,kl) -=
          (*this->singleSlater_->epsA())(i) + 
          (*this->singleSlater_->epsA())(j) - 2*this->rMu_; 
      }
    } // this->currentMat_ == CAB_PPTDA
    else if(this->currentMat_ == CBB_PPTDA){
    // ************************** //
    // ** BETA-BETA BLOCKS NYI ** //
    // ************************** //
    } // this->currentMat_ == CBB_PPTDA
    else if(this->currentMat_ == PPRPA_SINGLETS){
      auto iOff = this->nVAVA_LT_;
      auto iMetricScale = -1;
      // Place A Matrix
      for(auto a = 0, ab = 0; a < this->nVA_; a++      )
      for(auto b = 0        ; b <= a       ;  b++, ab++)
      for(auto c = 0, cd = 0; c < this->nVA_; c++      )
      for(auto d = 0        ; d <= c       ;  d++, cd++){
        double fact = 1.0;
        if(a == b) fact *= std::sqrt(0.5);
        if(c == d) fact *= std::sqrt(0.5);
        this->transDen_[iMat](ab,cd) = 
          fact * (this->mointegrals_->ABCD(a,c,b,d) +
                  this->mointegrals_->ABCD(a,d,b,c));  

        if(ab == cd) this->transDen_[iMat](ab,cd) +=
          (*this->singleSlater_->epsA())(a+this->nOA_) + 
          (*this->singleSlater_->epsA())(b+this->nOA_) - 2*this->rMu_; 
      }

      // Place C Matrix
      for(auto i = 0, ij = iOff; i < this->nOA_; i++      )
      for(auto j = 0           ; j <= i        ; j++, ij++)
      for(auto k = 0, kl = iOff; k < this->nOA_; k++      )
      for(auto l = 0           ; l <= k        ; l++, kl++){
        double fact = 1.0;
        if(i == j) fact *= std::sqrt(0.5);
        if(k == l) fact *= std::sqrt(0.5);
        this->transDen_[iMat](ij,kl) = 
          fact * (this->mointegrals_->IJKL(i,k,j,l) + 
                  this->mointegrals_->IJKL(i,l,j,k));

        if(ij == kl) this->transDen_[iMat](ij,kl) -=
          (*this->singleSlater_->epsA())(i) + 
          (*this->singleSlater_->epsA())(j) - 2*this->rMu_; 

        this->transDen_[iMat](ij,kl) *= iMetricScale;
      }

      // Place B Matricies
      for(auto a = 0, ab =    0; a < this->nVA_; a++      )
      for(auto b = 0           ; b <= a        ; b++, ab++)
      for(auto i = 0, ij = iOff; i < this->nOA_; i++      )
      for(auto j = 0           ; j <= i        ; j++, ij++) {
        double fact = 1.0;
        if(a == b) fact *= std::sqrt(0.5);
        if(i == j) fact *= std::sqrt(0.5);
        this->transDen_[iMat](ab,ij) = 
          fact * (this->mointegrals_->IAJB(i,a,j,b) + 
                  this->mointegrals_->IAJB(i,b,j,a));
        this->transDen_[iMat](ij,ab) = 
          iMetricScale*this->transDen_[iMat](ab,ij);
      }
    } // this->currentMat_ == PPRPA_SINGLETS
    else if(this->currentMat_ == A_PPTDA_SINGLETS){
      for(auto a = 0, ab = 0; a < this->nVA_; a++      )
      for(auto b = 0        ; b <= a       ;  b++, ab++)
      for(auto c = 0, cd = 0; c < this->nVA_; c++      )
      for(auto d = 0        ; d <= c       ;  d++, cd++){
        double fact = 1.0;
        if(a == b) fact *= std::sqrt(0.5);
        if(c == d) fact *= std::sqrt(0.5);
        this->transDen_[iMat](ab,cd) = 
          fact * (this->mointegrals_->ABCD(a,c,b,d) +
                  this->mointegrals_->ABCD(a,d,b,c));  
        if(ab == cd) this->transDen_[iMat](ab,cd) +=
          (*this->singleSlater_->epsA())(a+this->nOA_) + 
          (*this->singleSlater_->epsA())(b+this->nOA_) - 2*this->rMu_; 
      }
    } // this->currentMat_ == A_PPTDA_SINGLETS
    else if(this->currentMat_ == C_PPTDA_SINGLETS){
      for(auto i = 0, ij = 0; i < this->nOA_; i++      )
      for(auto j = 0        ; j <= i        ; j++, ij++)
      for(auto k = 0, kl = 0; k < this->nOA_; k++      )
      for(auto l = 0        ; l <= k        ; l++, kl++){
        double fact = 1.0;
        if(i == j) fact *= std::sqrt(0.5);
        if(k == l) fact *= std::sqrt(0.5);
        this->transDen_[iMat](ij,kl) = 
          fact * (this->mointegrals_->IJKL(i,k,j,l) + 
                  this->mointegrals_->IJKL(i,l,j,k));
        if(ij == kl) this->transDen_[iMat](ij,kl) -=
          (*this->singleSlater_->epsA())(i) + 
          (*this->singleSlater_->epsA())(j) - 2*this->rMu_; 
      }
    } // this->currentMat_ == C_PPTDA_SINGLETS
    else if(this->currentMat_ == PPRPA_TRIPLETS){
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
    } // this->currentMat_ == PPRPA_TRIPLETS
    else if(this->currentMat_ == A_PPTDA_TRIPLETS){
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
    } // this->currentMat_ == A_PPTDA_TRIPLETS
    else if(this->currentMat_ == C_PPTDA_TRIPLETS){
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
    } // this->currentMat_ == C_PPTDA_TRIPLETS

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
      msg += std::to_string(static_cast<int>(this->currentMat_));
      CErr(msg,this->fileio_->out);
    };

  //RealVecMap Eig(this->frequencies_[iMat].data(),this->nSingleDim_);
  //prettyPrint(this->fileio_->out,Eig,"Eig");
    // Cleanup LAPACK Memory
    delete[] WORK;
  }; // loop over iMat
}; // fullPPRPA

}; // namespace ChronusQ
