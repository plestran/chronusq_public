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
    typedef Eigen::Map<TMat> TCMMap;
    typedef Eigen::Map<TVec> TVecMap;

    // Dimension Variables
    size_t nSingleDim_;  ///< Dimension of the problem
    size_t nSek_;        ///< Number of desired solution vectors
    size_t nGuess_;      ///< Number of initial guess vectors
    size_t maxSubSpace_; ///< Maximum dimension of iterative subspace

    // Iteration related variables
    int maxMicroIter_;
    int maxMacroIter_;

    // Various tolerance values
    double residualTol_;

    // MetaData about QN Calculation
    int nMicroIter_;
    int nMacroIter_;
    int nTotalIter_;
   
    // In-Core storage of solution quantities
    TCMMatrix * solutionVecR_; ///< In-core storage of (right) solution vectors
    TCMMatrix * solutionVecL_; ///< In-core storage of (left)  solution vectors
    VectorXd  * omega_;        ///< Frequencies for solution vectors
    TCMMatrix * diag_;         ///< Diagonal elements of problem
    
  }; // class QNCallable

  template <typename T>
  class QuasiNewton2 {

    // Useful Eigen Typedefs
    typedef Eigen::Matrix<T,Dynamic,Dynamic,ColMajor> TMat;
    typedef Eigen::Matrix<T,Dynamic,1> TVec;
    typedef Eigen::Map<TMat> TCMMap;
    typedef Eigen::Map<TVec> TVecMap;


  public:
    QuasiNewton2(){;};

  }; // class QuasiNewton2
}; // namespace ChronusQ

#endif

