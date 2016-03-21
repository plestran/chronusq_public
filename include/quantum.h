/* 
 *  The Chronus Quantum (ChronusQ) software package is high-performace 
 *  computational chemistry software with a strong emphasis on explicitly 
 *  time-dependent and post-SCF quantum mechanical methods.
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
#ifndef INCLUDED_QUANTUM
#define INCLUDED_QUANTUM
#include <global.h>

namespace ChronusQ {

  enum DENSITY_TYPE {
    TOTAL,
    SPIN,
    ALPHA,
    BETA
  };
  /**
   *  Abstract Quantum Class. This class, in essence, takes all things from
   *  which one can write down a density matrix and generalizes the use of
   *  said density matrix in the computation of properties, etc.
   *
   *  As of now, Quantum only knows about densities expanded in a finite basis
   *  of atomic orbitals, but can in principle be generalized to form of a 
   *  density given the ability to contract with operators.
   */
  template<typename T>
  class Quantum {
    // Useful typedefs
    typedef Eigen::Matrix<T,Dynamic,Dynamic,ColMajor> TMatrix;

    protected:

    int nTCS_;
    // Pointers to TMatrix quantities that will store the density of
    // the quantum system in a finite basis
    std::unique_ptr<TMatrix> densityA_;
    std::unique_ptr<TMatrix> densityB_;


    template<typename Scalar, typename Left, typename Right>
    static Scalar OperatorTrace(const Left& A, const Right& B) {
      auto result = A.frobInner(B);
      return reinterpret_cast<Scalar(&)[2]>(result)[0];
    }
      


    template<typename Scalar, DENSITY_TYPE DenTyp, typename Op>
    Scalar OperatorSpinCombine(const Op& op) {
      double zero = 0.0;
      if(DenTyp == DENSITY_TYPE::TOTAL)
        if(!this->isClosedShell)
          return this->computePropertyAlpha<Scalar>(op) + 
                 this->computePropertyBeta<Scalar>(op);
        else
          return this->computePropertyAlpha<Scalar>(op);
      else if(DenTyp == DENSITY_TYPE::ALPHA)
        if(!this->isClosedShell)
          return this->computePropertyAlpha<Scalar>(op);
        else
          return 0.5*this->computePropertyAlpha<Scalar>(op);
      else if(DenTyp == DENSITY_TYPE::BETA)
        if(!this->isClosedShell)
          return this->computePropertyBeta<Scalar>(op);
        else
          return 0.5*this->computePropertyAlpha<Scalar>(op);
      else if(DenTyp == DENSITY_TYPE::SPIN)
        if(!this->isClosedShell)
          return this->computePropertyAlpha<Scalar>(op) -
                 this->computePropertyBeta<Scalar>(op);
        else
          return reinterpret_cast<Scalar(&)[2]>(zero)[0];
    }


    template<typename Scalar, typename Op> 
    Scalar computePropertyAlpha(const Op& op){
      return OperatorTrace<Scalar>((*this->densityA_),op);
    }
    template<typename Scalar, typename Op> 
    Scalar computePropertyBeta(const Op& op){
      return OperatorTrace<Scalar>((*this->densityB_),op);
    }
      
    public:
  
    bool isClosedShell;

    Quantum(){
      this->densityA_ = nullptr; 
      this->densityB_ = nullptr; 
      this->isClosedShell = false;
    };

    Quantum(unsigned int N) : Quantum(){
      this->allocDensity(N);
    };

    virtual void formDensity() = 0;
    void allocDensity(unsigned int N) {
      this->densityA_ = std::unique_ptr<TMatrix>(new TMatrix(N,N));
      if(!this->isClosedShell){
        this->densityB_ = std::unique_ptr<TMatrix>(new TMatrix(N,N));
      }
    };


    template<typename Scalar, DENSITY_TYPE DenTyp, typename Op> 
    Scalar computeProperty(const Op& op){
      return this->OperatorSpinCombine<Scalar,DenTyp>(op);
    }

    template<typename Scalar, DENSITY_TYPE DenTyp, typename Op> 
    std::vector<Scalar> computeProperty(const std::vector<Op>& op){
      std::vector<Scalar> results;
      for(typename std::vector<Op>::const_iterator it = op.begin(); 
          it != op.end(); ++it)
        results.push_back(this->computeProperty<Scalar,DenTyp,Op>(*it));
      return results;
    }


    inline int   nTCS(){ return this->nTCS_;};      
    inline TMatrix* densityA(){ return this->densityA_.get();};
    inline TMatrix* densityB(){ return this->densityB_.get();};
    // MPI Routines
    void mpiBCastDensity();
  };

  template<typename T>
  void Quantum<T>::mpiBCastDensity(){
  #ifdef CQ_ENABLE_MPI
    auto dataType = MPI_DOUBLE;
    if(typeid(T).hash_code() == typeid(dcomplex).hash_code())
      dataType = MPI_C_DOUBLE_COMPLEX;
  
    MPI_Bcast(this->densityA_->data(),this->densityA_->size(),dataType,0,
      MPI_COMM_WORLD);
    if(!this->isClosedShell && this->nTCS_ != 2)
      MPI_Bcast(this->densityB_->data(),this->densityB_->size(),dataType,0,
        MPI_COMM_WORLD);
  #endif
    ;
  };

}; // namespace ChronusQ



#endif
