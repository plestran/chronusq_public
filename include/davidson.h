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
#include <global.h>
#include <cerr.h>
#include <sdresponse.h>
#ifndef INCLUDED_DAVIDSON
#define INCLUDED_DAVIDSON
namespace ChronusQ {
/**
 *  A class to run Davidson Diagonalization written by David Williams-Young
 *
 */
  template <typename T>
  class Davidson {
    enum{
      CIS,
      RPA,
      CCSD
    };
    typedef Eigen::Matrix<T,Dynamic,Dynamic,RowMajor> TMat;
    typedef Eigen::Matrix<T,Dynamic,1> TVec;
    int     n_;          // Dimension of the problem (LDA)
    TMat*   mat_;        // The full matrix to be diagonalized (?)
    bool    hermetian_;  // Whether or not the problem is hemetian

    std::unique_ptr<TMat> guess_;      // Guess vectors
    std::unique_ptr<TVec> eigenvalues_;
    std::unique_ptr<TMat> eigenvector_;


    int     maxSubSpace_; // Maximum iterative subspace
    int     maxIter_;     // Maximum number of iterations (micro)
    int     MaxIter_;     // Maximum number of iterations (macro) 
    int     nSek_;        // Number of desired roots
    int     nGuess_;      // Number of given (or created) guesses

    void runMicro(ostream &output=cout);
    bool converged_;
    bool useLAPACK_;
    TMat (*AX_)(const TMat &, const TMat &) ;      // Function to form AX

    int method_;
    SDResponse * sdr_;


  public:
    // Run the Davidson
    inline void run(ostream &output=cout) {
      time_t currentTime;
      std::chrono::high_resolution_clock::time_point start,finish;
      std::chrono::duration<double> elapsed;
      output << bannerTop << endl << endl;
      time(&currentTime);
      output << "Davidson Diagonalization Started: " << ctime(&currentTime) << endl;
      this->printInfo(output);
      start = std::chrono::high_resolution_clock::now();
      this->runMicro(output);
      finish = std::chrono::high_resolution_clock::now();
      elapsed = finish-start;
      time(&currentTime);
      if(!this->converged_) CErr("Davidson Didn't Converge!",output);
      output << "Davidson Diagonalization Finished: " << ctime(&currentTime);
      output << "Time Elapsed: " << elapsed.count() << " sec" << endl;
      output << bannerEnd << endl << endl;;
    };

    // Default contructor (nothing interesting...)
    Davidson() {
      this->mat_         = nullptr;
      this->AX_          = NULL;
      this->guess_       = nullptr;
      this->eigenvalues_ = nullptr;
      this->eigenvector_ = nullptr;
      this->maxSubSpace_ = 250;
      this->maxIter_     = 20;
      this->MaxIter_     = 20;
      this->nSek_        = 0;
      this->nGuess_      = 0;
      this->n_           = 0;
      this->converged_   = false;
      this->hermetian_   = false;
      this->useLAPACK_   = true;
      this->method_      = -1;
      this->sdr_    = nullptr;
    }

    // Pass a ptr to a matrix to be diagonalized
    Davidson(TMat* A, int nSek) {
      this->maxSubSpace_ = 250;
      this->maxIter_     = 128;
      this->MaxIter_     = 20;
      this->converged_   = false;
      this->useLAPACK_   = false; // Use LAPACK by default
      this->mat_    = A;
      this->AX_     = NULL;
      this->nSek_   = nSek;
      this->nGuess_ = 2*nSek;
      this->n_      = A->cols();
      this->method_      = -1;
      this->sdr_    = nullptr;
      this->guess_  = 
        std::unique_ptr<TMat>(new TMat(this->n_,this->nGuess_));
      this->eigenvalues_ = 
        std::unique_ptr<TMat>(new TMat(this->nSek_,1));
      this->eigenvector_ = 
        std::unique_ptr<TMat>(new TMat(this->n_,this->nSek_));
/*
      this->guess_  = 
        std::make_shared<TMat>(this->n_,this->nGuess_);
      this->eigenvalues_ = 
        std::make_shared<TVec>(this->nSek_);
      this->eigenvector_ = 
        std::make_shared<TMat>(this->n_,this->nSek_);
*/
      this->hermetian_ = true; // Only supports Hermetian for time being

      (*this->guess_) = TMat::Identity(this->n_,this->nGuess_); // Identity guess (primitive)
    }

    // Pass a function to compute AX (dummy routine to test function passing
    Davidson(TMat (*AX)(const TMat&, const TMat&), TMat * A, int nSek, int N) {
      this->maxSubSpace_ = 250;
      this->maxIter_     = 128;
      this->MaxIter_     = 20;
      this->converged_   = false;
      this->useLAPACK_   = false; // Use LAPACK by default
      this->mat_    = A;
      this->AX_     = AX;
      this->nSek_   = nSek;
      this->nGuess_ = 2*nSek;
      this->n_      = N;
      this->method_      = -1;
      this->sdr_    = nullptr;
      this->guess_  = 
        std::unique_ptr<TMat>(new TMat(this->n_,this->nGuess_));
      this->eigenvalues_ = 
        std::unique_ptr<TMat>(new TMat(this->nSek_,1));
      this->eigenvector_ = 
        std::unique_ptr<TMat>(new TMat(this->n_,this->nSek_));
/*
      this->guess_  = 
        std::make_shared<TMat>(this->n_,this->nGuess_);
      this->eigenvalues_ = 
        std::make_shared<TVec>(this->nSek_);
      this->eigenvector_ = 
        std::make_shared<TMat>(this->n_,this->nSek_);
*/
      this->hermetian_ = true; // Only supports Hermetian for time being

      (*this->guess_) = TMat::Identity(this->n_,this->nGuess_); // Identity guess (primitive)
    }

    Davidson(SDResponse * SDR, int meth, int nSek){
      this->maxSubSpace_ = 250;
      this->maxIter_     = 128;
      this->MaxIter_     = 20;
      this->converged_   = false;
      this->useLAPACK_   = false; // Use LAPACK by default
      this->mat_    = nullptr;
      this->AX_     = NULL;
      this->nSek_   = nSek;
      this->nGuess_ = 2*nSek;
//    this->n_      = N;
      this->method_ = meth;
      this->sdr_    = SDR;
      this->guess_  = 
        std::unique_ptr<TMat>(new TMat(this->n_,this->nGuess_));
      this->eigenvalues_ = 
        std::unique_ptr<TMat>(new TMat(this->nSek_,1));
      this->eigenvector_ = 
        std::unique_ptr<TMat>(new TMat(this->n_,this->nSek_));
/*
      this->guess_  = 
        std::make_shared<TMat>(this->n_,this->nGuess_);
      this->eigenvalues_ = 
        std::make_shared<TVec>(this->nSek_);
      this->eigenvector_ = 
        std::make_shared<TMat>(this->n_,this->nSek_);
*/
      this->hermetian_ = true; // Only supports Hermetian for time being

      (*this->guess_) = TMat::Identity(this->n_,this->nGuess_); // Identity guess (primitive)

    }
    ~Davidson(){;};
    
    inline void printInfo(ostream &output=cout) {
      output << bannerTop << endl;
      output << "Davidson Diagonalization Settings:" << endl << endl;
 
      output << std::setw(50) << std::left << "  Dimension of the Matrix:" << this->n_ << endl;
      output << std::setw(50) << std::left << "  Number of Desired Roots:" << this->nSek_ << endl;
      output << std::setw(50) << std::left << "  Number of Initial Guess Vectors:" << this->nGuess_;
      if(this->nGuess_ == 2*this->nSek_) output << "    (default 2*NSek)";
      output << endl;
      output << std::setw(50) << std::left << "  Maximum Dimension of Iterative Subspace:"
           << this->maxSubSpace_ << endl;
      output << std::setw(50) << std::left << "  Maximum Number of Micro Iterations:"
           << this->maxIter_ << endl;
      output << std::setw(50) << std::left << "  Maximum Number of Macro Iterations:"
           << this->MaxIter_ << endl;
      output << std::setw(50) << std::left << "  Using an Hermetian algorithm?:";
      if(this->hermetian_) output << "Yes";
      else output << "No";
      output << endl;
      output << std::setw(50) << std::left << "  Using LAPACK to diagonalize subspace?:";
      if(this->useLAPACK_) output << "Yes";
      else output << "No";
      output << endl;
      output << std::setw(50) << std::left << "  Full Matrix Passed to for AX:";
      if(this->AX_==NULL) output << "Yes";
      else output << "No";
      output << endl;
 
      output << endl << bannerEnd << endl;
    }
    
  }; // class Davidson
} // namespace ChronusQ
#endif
