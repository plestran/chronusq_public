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
 * A class to setup and run various iterative Quasi-Newton type 
 * calculations such as:
 *   Davidson Diagonalization
 *   Quasi-Newton Linear Equation Solution
 */
template <typename T>
  class QuasiNewton{
    // Useful typedefs for Eigen templates
    typedef Eigen::Matrix<T,Dynamic,Dynamic,ColMajor> TMat;
    typedef Eigen::Matrix<T,Dynamic,1> TVec;

    bool isHermetian_;     // Hermetian Scheme
    bool doDiag_;          // Quasi-Newton Diagonalization (Davidson)
    bool doGEP_;           // Generalized Eigenproblem
    bool doLin_;           // Quasi-Newton Linear Equation Solve
    bool symmtrizedTrial_; // Symmetrized Trial Vectors (Kauczor et al. JCTC 7 (2010) p. 1610)
    bool debug_;           // Enables various debug options
    bool cleanupMem_;      // True if memory cleanup is required (if memory was allocated locally)
    bool isConverged_;     // True if the iterative calculation converged
    bool genGuess_;        // Generate Standard (identity) Guess

    int N_;                // Dimension of the problem
    int maxSubSpace_;      // Maximum dimension of the iterative subspace
    int maxIter_;          // Maximum number of miro iterations
    int MaxIter_;          // Maximum number of macro iterations
    int nSek_;             // Number of desired roots in the case of diagonalization (Davidson)
    int nGuess_;           // Numer of given (or to be generated) guess vectors
//  int method_;           // Post-SCF method being performed

    TMat * A_;                // Pointer to full matrix to be diagonalized
    TMat * solutionVector_;   // Solution vectors at current iteration (could be local or non-local)
    TVec * solutionValues_;   // Solution values (eigenvalues) at current iteration (could be local or non-local)
    std::unique_ptr<TMat> guess_;   // Guess vectors (always local copy, even if generated elsewhere)

    SDResponse * sdr_;     // Pointer to SDResponse object

    void runMicro(ostream &output=cout);
    void checkValid(ostream &output=cout);
    void loadDefaults();

    inline void allocGuess(){ 
      this->guess_ = std:unique_ptr<TMat>(new TMat(this->n_,this->nGuess_));
      if(this->genGuess_) this->identGuess();
    };
    inline void allocSolution(){
      this->solutionVector_ = new TMat(this->nSek_,1);
      this->solutionValues_ = new TMat(this->n_,this->nSek_);
      this->cleanupMem_ = true;
    };
    inline void identGuess(){
      (*this->guess_) = TMat::Identity(this->n_,this->nGuess_);
    };
    inline int stdSubSpace(){
      return std::min(6*this->nSek_,this->N_/2);
    };
    inline int stdNGuess(){
      return 2*this->nSek_;
    }

  public:
    /**
     *  Destructor
     *
     *  Only useful if memory for the solution was allocated as opposed
     *  to being passed in
     */ 
    ~QuasiNewton(){
      if(this->cleanupMem_){
        delete this->solutionVector_;
        delete this->solutionValues_;
      }
    }

    /**
     *  Default Constructor
     *
     *  Loads the default parameter values (see loadDefaults)
     */ 
    QuasiNewton(){ this->loadDefaults();};

    /**
     *  Constructor for debug comparison
     *  
     *  Full matrix is passed in so one can test any number of things:
     *    - Matrix-Vector product (sigma)
     *    - Quality of guess
     *    - Iteration comparison
     *
     *  This has no place in production code and should only be used for
     *  debugging. Storing the full matrix in memory completely defeats
     *  the purpose of iterative Quasi-Newton schemes.
     *
     *  Debug print is turned on by default
     */
    QuasiNewton(TMat* A, int nSek, ostream &output=cout) : QuasiNewton() {
      this->A_           = A;
      this->N_           = A.rows();
      this->nSek_        = nSek;
      this->nGuess_      = this->stdNGuess();
      this->maxSubSpace_ = this->stdSubSpace();
      this->debug_       = true;

      this->checkValid();
      this->allocGuess();
      this->allocSolution();
    };

    /**
     *  Constructor for Quasi-Newton based on a SDResponse Object
     *
     *  All of the parameters are set based on the SDResponse object
     *  that is passed in. Currently, QuasiNewton knows about the
     *  following from SDResponse:
     *   -  Configuration Interaction Singles (CIS)
     *   -  Random Phase Approximation        (RPA)
     *   
     *  Currently, the following developments are under-way and
     *  QuasiNewton will know how to handle them soon enough
     *   -  Standard Frequency Dependent Linear Response (SFDLR)
     *   -  Damped Frequency Dependent Linear Response   (DFDLR)
     */ 
    QuasiNewton(SDResponse * SDR) : QuasiNewton(){
      this->sdr_              = SDR; 
      this->N_                = SDR->nSingleDim();
      this->nSek_             = SDR->nSek();
      this->nGuess_           = SDR->nGuess();
      this->solutionValues_   = SDR->omega();
      this->solutionVector_   = SDR->transDen();
      this->symmetrizedTrial_ = (SDR->iMeth() == SDResponse::RPA);
      this->doGEP_            = (SDR->iMeth() == SDResponse::RPA);
      this->genGuess_         = (this->nGuess_ == 0);

      if(this->genGuess_) this->nGuess_ = this->stdNGuess();
      this->checkValid();
      this->allocGuess();
      this->allocSolution();
      if(!this->genGuess_) *this->guess_ = *SDR->davGuess();
    };

   inline TVec* eigenValues(){return this->solutionValues_;};
   inline TMat* eigenVector(){return this->solutionVector_;};
   inline void run(ostream &output=cout){
     time_t currentTime;
     std::chrono::high_resolution_clock::time_pint start,finish;
     std::chrono::duration<double> elapsed;
     output << bannerTop << endl << endl;
     time(&currentTime);
     outout << "Quasi-Newton Calculation Started: " << ctime(&currentTime) << endl;
     this->printInfo(output);
     start = std::chrono::high_resolution_clock::now();
     this->runMicro(output);
     finish = std::chrono::high_resolution_clock::now();
     elapsed = finish - start;
     time(&currentTime);
     if(!this->isConverged_) CErr("Quasi-Newton Failed to Converge Within Given Criteria!",output);
     output << "Quasi-Newton Calculation Finished: " << ctime(&currentTime);
     output << "Time Elapsed: " << elapsed.count() << " sec" << endl;
     output << bannerEnd << endl << endl;
   }
  
  } // class QuasiNewton

  void QuasiNewton::loadDefaults(){
    this->isHermetian_      = true;  // Default to Hermetian Scheme
    this->doDiag_           = true;  // Defualt to diagonalization
    this->doGEP_            = false; // Default standard eigenproblem if doing a diagonalization
    this->doLin_            = false; // Mildly redundent, as it will currently always be the opposite of doDiag
    this->symmetrizedTrial_ = false; // Default to unmodified trial vectors
    this->debug_            = false; // Default to terse output
    this->cleanupMem_       = false; // Assume that the space for results was allocated elsewhere
    this->isConverged_      = false; // Obviously the calculation doesn't begin converged...
    this->genGuess_         = true;  // Defualt generate identity guess

    // Zero out all of the size dependent quantities
    this->N_                = 0;
    this->maxSubSpace_      = 0;
    this->nSek_             = 0;
    this->nGuess_           = 0; 

    // Defaults for # of interations
    this->maxIter_          = 128;
    this->MaxIter_          = 20;

    // Initialize all pointers to some varient of NULL
    this->A_                = NULL;
    this->solutionVector_   = NULL;
    this->solutionValues_   = NULL;
    this->guess_            = nullptr;
    this->sdr_              = NULL;
  }

  void QuasiNewton::checkValid(ostream &output){
    if(this->A_ != NULL){
      if(this->A_->rows() != this->A_->cols())
        CErr("Quasi-Newton only supported for square problems!")
    }
    if(this->nGuess_ >= this->maxSubSpace_)
      CErr("Number of initial guess vectors exceeds maximum dimension of iterative subspace");
  }
} // namespace ChronusQ
