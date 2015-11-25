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
    char JOBZ = 'V';
    char UPLO = 'U';
    int INFO;
    int LWORK = 6*this->nSingleDim_;
    double * WORK  = new double[LWORK]; 
    double * W     = new double[this->nSingleDim_];

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

    int N = this->nSingleDim_;
    if(this->doTDA_)
      dsyev_(&JOBZ,&UPLO,&N,this->transDen_[iMat].data(),&N,
        this->frequencies_[iMat].data(),WORK,&LWORK,&INFO);
    else {
//    prettyPrint(this->fileio_->out,this->transDen_[iMat],"Mat");
      double * COPY = new double[this->nSingleDim_*this->nSingleDim_]; 
      std::memcpy(COPY,this->transDen_[iMat].data(),this->nSingleDim_*this->nSingleDim_*sizeof(double));
      double * WI = new double[this->nSingleDim_];
      std::memset(WI,0.0,this->nSingleDim_*sizeof(double));
      char JOBVL = 'N';

      dgeev_(&JOBVL,&JOBZ,&N,COPY,&N,this->frequencies_[iMat].data(),WI,
        this->transDen_[iMat].data(),&N,this->transDen_[iMat].data(),&N,
        WORK,&LWORK,&INFO);

      std::sort(this->frequencies_[iMat].data(),
                this->frequencies_[iMat].data()+
                this->frequencies_[iMat].size());
    }
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
    delete[] W;

  } 

}; // fullFOPPA (T = double)
}; // namespace ChronusQ
