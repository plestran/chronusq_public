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

    int  nTCS_;
    int  maxMultipole_;
    bool isScattered_;


    // Pointers to TMatrix quantities that will store the density of
    // the quantum system in a finite basis

    // The Alpha and Beta Components of the Density matrix. In the case
    // of spinors, onePDMA holds the entire density (AA,AB,BA,BB)
    std::unique_ptr<TMatrix> onePDMA_;
    std::unique_ptr<TMatrix> onePDMB_;

    // The Density in the basis of the identity and the Pauli matricies 
    std::unique_ptr<TMatrix> onePDMScalar_;
    std::unique_ptr<TMatrix> onePDMMz_;
    std::unique_ptr<TMatrix> onePDMMy_;
    std::unique_ptr<TMatrix> onePDMMx_;

    std::array<double,3> elecDipole_;
    std::array<std::array<double,3>,3> elecQuadpole_;
    std::array<std::array<double,3>,3> elecTracelessQuadpole_;
    std::array<std::array<std::array<double,3>,3>,3> elecOctpole_;


    template<typename Scalar, typename Left, typename Right>
    static Scalar OperatorTrace(const Left& A, const Right& B) {
      auto result = A.frobInner(B);
      return reinterpret_cast<Scalar(&)[2]>(result)[0];
    }
      


    template<typename Scalar, DENSITY_TYPE DenTyp, typename Op>
    Scalar OperatorSpinCombine(const Op& op) {
      double zero = 0.0;
      if(DenTyp == DENSITY_TYPE::TOTAL){
        if(!this->isClosedShell)
          return this->computePropertyAlpha<Scalar>(op) + 
                 this->computePropertyBeta<Scalar>(op);
        else
          return this->computePropertyAlpha<Scalar>(op);
      } else if(DenTyp == DENSITY_TYPE::ALPHA) {
        if(!this->isClosedShell)
          return this->computePropertyAlpha<Scalar>(op);
        else
          return 0.5*this->computePropertyAlpha<Scalar>(op);
      } else if(DenTyp == DENSITY_TYPE::BETA) {
        if(!this->isClosedShell)
          return this->computePropertyBeta<Scalar>(op);
        else
          return 0.5*this->computePropertyAlpha<Scalar>(op);
      } else if(DenTyp == DENSITY_TYPE::SPIN) {
        if(!this->isClosedShell)
          return this->computePropertyAlpha<Scalar>(op) -
                 this->computePropertyBeta<Scalar>(op);
        else
          return reinterpret_cast<Scalar(&)[2]>(zero)[0];
      }
    }


    template<typename Scalar, typename Op> 
    Scalar computePropertyAlpha(const Op& op){
      return OperatorTrace<Scalar>((*this->onePDMA_),op);
    }
    template<typename Scalar, typename Op> 
    Scalar computePropertyBeta(const Op& op){
      return OperatorTrace<Scalar>((*this->onePDMB_),op);
    }
      
    public:
  
    bool isClosedShell;

    Quantum(){
      this->onePDMA_ = nullptr; 
      this->onePDMB_ = nullptr; 
      this->isClosedShell = false;
      this->nTCS_ = 1;
      this->maxMultipole_ = 3;
      this->isScattered_ = false;


      this->clearElecMultipole();
    };

    Quantum(unsigned int N) : Quantum(){
      this->allocDensity(N);
    };

    Quantum(const Quantum &other) :  // Copy Constructor
      elecDipole_(other.elecDipole_),
      elecQuadpole_(other.elecQuadpole_),
      elecTracelessQuadpole_(other.elecTracelessQuadpole_),
      elecOctpole_(other.elecOctpole_),
      nTCS_(other.nTCS_),
      isClosedShell(other.isClosedShell),
      maxMultipole_(other.maxMultipole_) {

      this->onePDMA_ = std::unique_ptr<TMatrix>(
          new TMatrix(*other.onePDMA_)
        );
      if(!this->isClosedShell && this->nTCS_ != 2)
        this->onePDMB_ = std::unique_ptr<TMatrix>(
            new TMatrix(*other.onePDMB_)
          );
    }

    template<typename U>
    Quantum(const U&);

    virtual void formDensity() = 0;
    inline void allocDensity(unsigned int N) {
      if(!this->isScattered_) {
        this->onePDMA_ = 
          std::unique_ptr<TMatrix>(
              new TMatrix(this->nTCS_ * N,this->nTCS_ * N)
          );
        if(!this->isClosedShell){
          this->onePDMB_ = 
            std::unique_ptr<TMatrix>(
                new TMatrix(this->nTCS_ * N,this->nTCS_ * N)
            );
        }
      } else {
        this->onePDMScalar_ = std::unique_ptr<TMatrix>(new TMatrix(N,N));
        if((this->nTCS_ == 1 && !this->isClosedShell) || this->nTCS_ == 2)
          this->onePDMMz_ = std::unique_ptr<TMatrix>(new TMatrix(N,N));
        if(this->nTCS_ == 2) {
          this->onePDMMx_ = std::unique_ptr<TMatrix>(new TMatrix(N,N));
          this->onePDMMy_ = std::unique_ptr<TMatrix>(new TMatrix(N,N));
        }
      }
    };

    inline void alloc(unsigned int N){
      this->allocDensity(N);
    };

    inline void clearElecMultipole(){
      for(auto iXYZ = 0; iXYZ < 3; iXYZ++){
        this->elecDipole_[iXYZ] = 0.0; 
        for(auto jXYZ = 0; jXYZ < 3; jXYZ++){
          this->elecQuadpole_[iXYZ][jXYZ] = 0.0; 
          this->elecTracelessQuadpole_[iXYZ][jXYZ] = 0.0; 
          for(auto kXYZ = 0; kXYZ < 3; kXYZ++)
            this->elecOctpole_[iXYZ][jXYZ][kXYZ] = 0.0; 
        }
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
    #include <quantum/quantum_stdproperties.h>


    void scatterDensity();
    void gatherDensity();
    void complexMyScale();


    inline void setMaxMultipole(int i){ this->maxMultipole_ = i;   };
    inline void setNTCS(int i){         this->nTCS_ = i;           };

    inline int   nTCS(){ return nTCS_;};      
    inline int maxMultipole(){ return maxMultipole_;};
    inline TMatrix* onePDMA(){ return onePDMA_.get();};
    inline TMatrix* onePDMB(){ return onePDMB_.get();};
    inline TMatrix* densityA(){ return onePDMA_.get();};
    inline TMatrix* densityB(){ return onePDMB_.get();};

    inline std::array<double,3> elecDipole(){ return elecDipole_; };
    inline std::array<std::array<double,3>,3> elecQuadpole(){ 
      return elecQuadpole_; 
    };
    inline std::array<std::array<double,3>,3> elecTracelessQuadpole(){ 
      return elecTracelessQuadpole_; 
    };
    inline std::array<std::array<std::array<double,3>,3>,3> elecOctpole(){ 
      return elecOctpole_; 
    };


    // MPI Routines
    void mpiBCastDensity();
  };

  template<typename T>
  void Quantum<T>::mpiBCastDensity(){
  #ifdef CQ_ENABLE_MPI
    auto dataType = MPI_DOUBLE;
    if(typeid(T).hash_code() == typeid(dcomplex).hash_code())
      dataType = MPI_C_DOUBLE_COMPLEX;
  
    MPI_Bcast(this->onePDMA_->data(),this->onePDMA_->size(),dataType,0,
      MPI_COMM_WORLD);
    if(!this->isClosedShell && this->nTCS_ != 2)
      MPI_Bcast(this->onePDMB_->data(),this->onePDMB_->size(),dataType,0,
        MPI_COMM_WORLD);
  #endif
    ;
  };

  template<typename T>
  void Quantum<T>::scatterDensity(){
    if(this->isScattered_) return;

    this->isScattered_ = true;
    size_t currentDim = this->onePDMA_->cols();
    // Allocate new scattered densities
    this->allocDensity(currentDim / this->nTCS_);

    if(this->nTCS_ == 1 && !this->isClosedShell) {

      (*this->onePDMScalar_) = (*this->onePDMA_) + (*this->onePDMB_);
      (*this->onePDMMz_)     = (*this->onePDMA_) - (*this->onePDMB_);

      // Deallocate space
      this->onePDMA_.reset();
      this->onePDMB_.reset();

    } else if(this->nTCS_ == 2) {

      for(auto I = 0, i = 0; I < currentDim; I += this->nTCS_, i++)
      for(auto J = 0, j = 0; J < currentDim; J += this->nTCS_, j++){

        (*this->onePDMScalar_)(i,j) = 
          (*this->onePDMA_)(I,J) + (*this->onePDMA_)(I+1,J+1);
        (*this->onePDMMz_)(i,j) = 
          (*this->onePDMA_)(I,J) - (*this->onePDMA_)(I+1,J+1);
        (*this->onePDMMx_)(i,j) = 
          (*this->onePDMA_)(I+1,J) + (*this->onePDMA_)(I,J+1);

        // In the case of Real Orbitals, My is a pure imaginary
        // Hermetian matrix stored as an anti-symmetric 
        // real matrix
        (*this->onePDMMy_)(i,j) = 
          (*this->onePDMA_)(I,J+1) - (*this->onePDMA_)(I+1,J);

      };
      /*
      Eigen::Map<TMatrix,0,Eigen::Stride<Dynamic,Dynamic> > 
        PAA(this->onePDMA_->data(),
            currentDim/this->nTCS_, currentDim/this->nTCS_,
            Eigen::Stride<Dynamic,Dynamic>(2,2));
      Eigen::Map<TMatrix,0,Eigen::Stride<Dynamic,Dynamic> > 
        PBA(this->onePDMA_->data() + 1,
            currentDim/this->nTCS_, currentDim/this->nTCS_,
            Eigen::Stride<Dynamic,Dynamic>(2,2));
      Eigen::Map<TMatrix,0,Eigen::Stride<Dynamic,Dynamic> > 
        PAB(this->onePDMA_->data() + currentDim/this->nTCS_,
            currentDim/this->nTCS_, currentDim/this->nTCS_,
            Eigen::Stride<Dynamic,Dynamic>(2,2));
      Eigen::Map<TMatrix,0,Eigen::Stride<Dynamic,Dynamic> > 
        PBB(this->onePDMA_->data() + (currentDim/this->nTCS_) + 1,
            currentDim/this->nTCS_, currentDim/this->nTCS_,
            Eigen::Stride<Dynamic,Dynamic>(2,2));

      this->onePDMScalar_->noalias() = PAA + PBB;
      this->onePDMMz_->noalias()     = PAA - PBB;
      this->onePDMMy_->noalias()     = PBA - PAB;
      this->onePDMMx_->noalias()     = PBA + PAB;
      */

      this->complexMyScale();

      this->onePDMA_.reset();

    };
  };

  template<typename T>
  void Quantum<T>::gatherDensity(){
    if(!this->isScattered_) return;

    this->isScattered_ = false;
    size_t currentDim = this->onePDMScalar_->cols();
    // Allocate new scattered densities
    this->allocDensity(currentDim);

    if(this->nTCS_ == 1 && !this->isClosedShell) {

      this->onePDMA_->noalias() = (*this->onePDMScalar_) + (*this->onePDMMz_);
      this->onePDMB_->noalias() = (*this->onePDMScalar_) - (*this->onePDMMz_);

      (*this->onePDMA_) *= 0.5;
      (*this->onePDMB_) *= 0.5;

      // Deallocate space
      this->onePDMScalar_.reset();
      this->onePDMMz_.reset();

    } else if(this->nTCS_ == 2) {

      // Since 
      //   My = i(PBA - PAB)
      //   PAB = Mx - i * My
      //   PBA = Mx + i * My
      // "My" is scaled by "i" before entering the reconstruction
      //
      // ** Note that this scaling is a dummy call for double 
      // precision objects and the sign is accounted for implicitly
      // through a flip in sign in the reconstruction **
      this->complexMyScale();

      for(auto I = 0, i = 0; i < currentDim; I += this->nTCS_, i++)
      for(auto J = 0, j = 0; j < currentDim; J += this->nTCS_, j++){

        (*this->onePDMA_)(I,J) = 
          (*this->onePDMScalar_)(i,j) + (*this->onePDMMz_)(i,j);
        (*this->onePDMA_)(I+1,J+1) = 
          (*this->onePDMScalar_)(i,j) - (*this->onePDMMz_)(i,j);

        if(typeid(T).hash_code() == typeid(dcomplex).hash_code()){
          (*this->onePDMA_)(I,J+1) = 
            (*this->onePDMMx_)(i,j) - (*this->onePDMMy_)(i,j); 
          (*this->onePDMA_)(I+1,J) = 
            (*this->onePDMMx_)(i,j) + (*this->onePDMMy_)(i,j); 
        } else {
          // Sign flip viz complex case because there is an implied
          // "i" infront of the pure imaginary My
          (*this->onePDMA_)(I,J+1) = 
            (*this->onePDMMx_)(i,j) + (*this->onePDMMy_)(i,j); 
          (*this->onePDMA_)(I+1,J) = 
            (*this->onePDMMx_)(i,j) - (*this->onePDMMy_)(i,j); 
        }
      };

      (*this->onePDMA_) *= 0.5;

      // Deallocate Space
      this->onePDMScalar_.reset();
      this->onePDMMz_.reset();
      this->onePDMMy_.reset();
      this->onePDMMx_.reset();

    };
  };



}; // namespace ChronusQ



#endif
