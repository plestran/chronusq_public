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
#ifndef INCLUDED_QUASINEWTON2 
#define INCLUDED_QUATINEWTON2
#include <global.h>
#include <cerr.h>

namespace ChronusQ { 
  template<typename T>
  class QNCallable {
    // Useful Eigen Typedefs
    typedef Eigen::Matrix<T,Dynamic,Dynamic,ColMajor> TMat;
    typedef Eigen::Matrix<T,Dynamic,1> TVec;
    typedef Eigen::Map<TMat> TMap;
    typedef Eigen::Map<TVec> TVecMap;

    // Dimension Variables
    size_t nSingleDim_;  ///< Dimension of the problem
    size_t nSek_;        ///< Number of desired solution vectors
    size_t nGuess_;      ///< Number of initial guess vectors

    // Booleans for QN setup
    bool generateGuess_;  ///< (?) Generate identity guess
    bool allocSolution_;  ///< (?) Allocate space for solution
    bool needsLeft_;      ///< (?) Left vectors need to be taken into account
    bool isFullMatrix_;   ///< (?) Full matrix passes (DEBUG OPTION)
    bool directDiag_;     ///< (?) Arithmetic with direct access of diagonals
   
    // In-Core storage
    TMat        * fullMatrix_;   ///< Full Matrix (DEBUG PURPOSES)
    TMat        * fullMetric_;   ///< Full Metric (DEBUG PURPOSES)
    TMat        * solutionVecR_; ///< (right) solution vectors
    TMat        * solutionVecL_; ///< (left)  solution vectors
    VectorXd    * omega_;        ///< Frequencies for solution vectors
    TMat        * diag_;         ///< Diagonal elements of problem

    // Disk storage
    H5::H5File  * scratchFile_; ///< Scratch file
    H5::DataSet * guessFile_;   ///< File to store guess vectors

  public:
    /**
     *  Defualt Constructor for QNCallable
     *
     *  Default constructor loads defaults
     */ 
    QNCallable() {
      this->nSingleDim_   = 0;
      this->nSek_         = 0;
      this->nGuess_       = 0;
      this->fullMatrix_   = NULL;
      this->solutionVecR_ = NULL;
      this->solutionVecL_ = NULL;
      this->omega_        = NULL;
      this->diag_         = NULL;
      this->scratchFile_  = NULL;
      this->guessFile_    = NULL;

      this->generateGuess_ = false;
      this->allocSolution_ = false;
      this->needsLeft_     = false;
      this->isFullMatrix_  = false;
      this->directDiag_    = false;
    };

    QNCallable(TMat *A, int NSek, int NGuess) : QNCallable(){
      this->fullMatrix_ = A;
      this->nSek_       = NSek;
      this->nGuess_     = NGuess;

      this->nSingleDim_    = A->rows();
      this->generateGuess_ = true;
      this->allocSolution_ = true;
      this->isFullMatrix_  = true;
    };

    QNCallable(TMat *A, int NSek) : QNCallable(A,NSek,0) {
      this->nGuess_ = 2*NSek;
    };

    QNCallable(TMat *A, int NSek, int NGuess, H5::H5File *scr) :
      QNCallable(A,NSek,NGuess) {
      this->scratchFile_ = scr;
    }

    inline void initQN() {
      if(this->generateGuess_) this->generateGuess();
      if(this->allocSolution_) this->allocSolution();
    }

    /** Function Declarations **/

    // Setters
    inline void setNSek(int n)  { this->nSek_   = n;};
    inline void setNGuess(int n){ this->nGuess_ = n;};

    // Getters
    inline int nGuess()     { return this->nGuess_    ; };
    inline int nSek()       { return this->nSek_      ; };
    inline int nSingleDim() { return this->nSingleDim_; };

    inline VectorXd * omega()        { return this->omega_       ; };
    inline TMat     * solutionVecR() { return this->solutionVecR_; };
    inline TMat     * solutionVecL() { return this->solutionVecL_; };


    virtual void linearTrans(TMap &,TMap &,TMap &,TMap &,TMap &,TMap &) = 0;
    virtual void formGuess() = 0;
    virtual void formDiag()  = 0;

    
  }; // class QNCallable

  // Enums for Job Control
  enum QNProblemType {
    DIAGONALIZATION,
    LINEAR_SOLVE
  };

  enum QNMatrixType {
    HERMETIAN,
    HERMETIAN_GEP,
    NON_HERMETIAN
  };

  enum QNGuessType {
    RESIDUAL_DAVIDSON,
    RESIDUAL_OLSEN
  };

  enum QNSpecialAlgorithm {
    NOT_SPECIAL,
    SYMMETRIZED_TRIAL
  };

  template <typename T>
  class QuasiNewton2 {

    // Useful Eigen Typedefs
    typedef Eigen::Matrix<T,Dynamic,Dynamic,ColMajor> TMat;
    typedef Eigen::Matrix<T,Dynamic,1> TVec;
    typedef Eigen::Map<TMat> TCMMap;
    typedef Eigen::Map<TVec> TVecMap;

    QNCallable<T> * qnObj_;

    // Dimension Variables
    size_t maxSubSpace_; ///< Maximum dimension of iterative subspace

    // Iteration related variables
    int maxMicroIter_;
    int maxMacroIter_;

    // Various tolerance values
    double residualTol_;

    // Job Control
    QNProblemType      problemType_;
    QNMatrixType       matrixType_; 
    QNGuessType        guessType_;
    QNSpecialAlgorithm specialAlgorithm_;

    // MetaData about QN Calculation
    int  nMicroIter_;
    int  nMacroIter_;
    int  nTotalIter_;
    bool isConverged_;


  public:
    #include <qn_constructors.h>

  }; // class QuasiNewton2
}; // namespace ChronusQ

#endif

