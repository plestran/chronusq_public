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
#include <tools.h>
#include <memory.h>

namespace ChronusQ {

  enum DENSITY_TYPE {
    TOTAL,
    SPIN,
    ALPHA,
    BETA,
    MZ,MX,MY
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
    typedef Eigen::Map<TMatrix> TMap;

    protected:

    int  nTCS_;
    int  maxMultipole_;
//  bool isScattered_;
    bool isAllocated_;

    CQMemManager * memManager_;


    // Pointers to TMatrix quantities that will store the density of
    // the quantum system in a finite basis

    // The Alpha and Beta Components of the Density matrix. In the case
    // of spinors, onePDMA holds the entire density (AA,AB,BA,BB)
//  std::unique_ptr<TMap> onePDMA_;
//  std::unique_ptr<TMap> onePDMB_;

    // The Density in the basis of the identity and the Pauli matricies 
    std::unique_ptr<TMap> onePDMScalar_;
    std::unique_ptr<TMap> onePDMMz_;
    std::unique_ptr<TMap> onePDMMy_;
    std::unique_ptr<TMap> onePDMMx_;

    std::vector<TMap*> onePDM_;
    
    /*
    // Breaks Up Scattered Density into Re/Im parts
    std::unique_ptr<RealMatrix> ReOnePDMScalar_;
    std::unique_ptr<RealMatrix> ReOnePDMMx_;
    std::unique_ptr<RealMatrix> ReOnePDMMy_;
    std::unique_ptr<RealMatrix> ReOnePDMMz_;
    std::unique_ptr<RealMatrix> ImOnePDMScalar_;
    std::unique_ptr<RealMatrix> ImOnePDMMx_;
    std::unique_ptr<RealMatrix> ImOnePDMMy_;
    std::unique_ptr<RealMatrix> ImOnePDMMz_;
    */


    std::array<double,3> elecDipole_;
    std::array<std::array<double,3>,3> elecQuadpole_;
    std::array<std::array<double,3>,3> elecTracelessQuadpole_;
    std::array<std::array<std::array<double,3>,3>,3> elecOctpole_;
    double Sx_, Sy_, Sz_, Ssq_;


    template<typename Scalar, typename Left, typename Right>
    static Scalar OperatorTrace(const Left& A, const Right& B) {
      auto result = A.frobInner(B);
      return reinterpret_cast<Scalar(&)[2]>(result)[0];
    }
      


    /*
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
    */
    template<typename Scalar, DENSITY_TYPE DenTyp, typename Op>
    Scalar OperatorSpinCombine(const Op& op) {
      double zero = 0.0;
      bool isReal = typeid(T).hash_code() == typeid(double).hash_code();

/*
      if(this->nTCS_ == 1 && this->isClosedShell){
        if(DenTyp == DENSITY_TYPE::TOTAL)
          return OperatorTrace<Scalar>((*this->onePDMA_),op);
        else
          return reinterpret_cast<Scalar(&)[2]>(zero)[0];
      } else {
*/
        if(DenTyp == DENSITY_TYPE::TOTAL)
          return OperatorTrace<Scalar>((*this->onePDMScalar_),op);
        else if(DenTyp == DENSITY_TYPE::SPIN || DenTyp == DENSITY_TYPE::MZ)
          return OperatorTrace<Scalar>((*this->onePDMMz_),op);
        else if(this->nTCS_ == 1)
          return reinterpret_cast<Scalar(&)[2]>(zero)[0];
        else {
          if(DenTyp == DENSITY_TYPE::MX)
            return OperatorTrace<Scalar>((*this->onePDMMx_),op);
          else if(isReal)
            return reinterpret_cast<Scalar(&)[2]>(zero)[0];
          else
            return OperatorTrace<Scalar>((*this->onePDMMy_),op);
        }
//    }
    }


/*
    template<typename Scalar, typename Op> 
    Scalar computePropertyAlpha(const Op& op){
      return OperatorTrace<Scalar>((*this->onePDMA_),op);
    }
    template<typename Scalar, typename Op> 
    Scalar computePropertyBeta(const Op& op){
      return OperatorTrace<Scalar>((*this->onePDMB_),op);
    }
*/
      
    public:
  
    bool isClosedShell;

    Quantum(){
    //this->onePDMA_ = nullptr; 
    //this->onePDMB_ = nullptr; 
      this->isClosedShell = false;
      this->nTCS_ = 1;
      this->maxMultipole_ = 3;
//    this->isScattered_ = false;
      this->isAllocated_ = false;

      cout << "In Quantum Constructor" << endl;
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
      maxMultipole_(other.maxMultipole_),
//    isScattered_(other.isScattered_),
      memManager_(other.memManager_) {

      this->isAllocated_ = false;
      auto NB = other.onePDMScalar_->rows(); 
//    auto NBT = other.onePDMScalar_->rows(); 
//    auto NBTSq = NBT*NBT;
//    this->alloc(NBT/this->nTCS_);
      this->alloc(NB);
     
     
/*
//    (*this->onePDMA_) = (*other.onePDMA_);
     
      if(!this->isClosedShell || this->nTCS_ == 2){
        (*this->onePDMScalar_) = (*other.onePDMScalar_);
        (*this->onePDMMz_)     = (*other.onePDMMz_);
        if(this->nTCS_ == 1)
          (*this->onePDMB_) = (*other.onePDMB_);
        else {
          (*this->onePDMMx_) = (*other.onePDMMx_);
          (*this->onePDMMy_) = (*other.onePDMMy_);
        } // NTCS check
      } // if not RHF/KS
*/
      for(auto IDen = 0; IDen < this->onePDM_.size(); IDen++)
        *this->onePDM_[IDen] = *other.onePDM_[IDen];
    }

    template<typename U>
    Quantum(const U&);

    // Link up to all of the other worker classes
    inline void communicate(CQMemManager &memManager){
      cout << "In Quantum Communicate" << endl;
      this->memManager_ = &memManager;
    }

    virtual void formDensity() = 0;
    inline void allocDensity(unsigned int N) {
      auto NSq = N*N;

/*
      this->onePDMA_ = 
        std::unique_ptr<TMap>(
            new TMap(this->memManager_->template malloc<T>(NBSq),NB,NB)
        );
      this->onePDMA_->setZero();
      if(!this->isClosedShell && this->nTCS_ == 1){
        this->onePDMB_ = 
          std::unique_ptr<TMap>(
            new TMap(this->memManager_->template malloc<T>(NBSq),NB,NB)
          );
        this->onePDMB_->setZero();
      }
*/

      this->onePDMScalar_ = std::unique_ptr<TMap>(
          new TMap(this->memManager_->template malloc<T>(NSq),N,N));

      this->onePDMScalar_->setZero();
      this->onePDM_.emplace_back(this->onePDMScalar_.get());

/*
      if(this->nTCS_ == 1 and this->isClosedShell) {
        this->onePDM_.emplace_back(this->onePDMA_.get());
      } else {
        this->onePDM_.emplace_back(this->onePDMScalar_.get());
      }
*/

      if(!this->isClosedShell || this->nTCS_ == 2){
        this->onePDMMz_ = std::unique_ptr<TMap>(
            new TMap(this->memManager_->template malloc<T>(NSq),N,N));

        this->onePDMMz_->setZero();
        this->onePDM_.emplace_back(this->onePDMMz_.get());
      }
      if(this->nTCS_ == 2) {
        this->onePDMMx_ = std::unique_ptr<TMap>(
            new TMap(this->memManager_->template malloc<T>(NSq),N,N));
        this->onePDMMy_ = std::unique_ptr<TMap>(
            new TMap(this->memManager_->template malloc<T>(NSq),N,N));

        this->onePDMMy_->setZero();
        this->onePDMMx_->setZero();
        this->onePDM_.emplace_back(this->onePDMMy_.get());
        this->onePDM_.emplace_back(this->onePDMMx_.get());
      }
    };

    inline void alloc(unsigned int N){
      if(this->isAllocated_) return;
      this->allocDensity(N);
      this->isAllocated_ = true;
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


    template<typename Op> 
    static void spinScatter( Op &, std::vector<std::reference_wrapper<Op>> &);
    template<typename Op> 
    static void spinScatter( Op &, Op &, 
        std::vector<std::reference_wrapper<Op>> &);

    template<typename Op> 
    static void spinGather( Op & , std::vector<std::reference_wrapper<Op>> &);
    template<typename Op> 
    static void spinGather( Op & , Op &, 
        std::vector<std::reference_wrapper<Op>> &);

//  void scatterDensity();
//  void gatherDensity();

//  template<typename Op>
//  static void complexMyScale(Op &);

//  static T ComplexScale();

    inline void setMaxMultipole(int i){ this->maxMultipole_ = i;   };
    inline void setNTCS(int i){         this->nTCS_ = i;           };

    inline int   nTCS(){ return nTCS_;};      
    inline int maxMultipole(){ return maxMultipole_;};
//  inline TMap* onePDMA(){ return onePDMA_.get();};
//  inline TMap* onePDMB(){ return onePDMB_.get();};
//  inline TMap* densityA(){ return onePDMA_.get();};
//  inline TMap* densityB(){ return onePDMB_.get();};
  
    inline std::vector<TMap*>& onePDM(){ return this->onePDM_;};
    inline TMap* onePDMScalar(){ return onePDMScalar_.get();};
    inline TMap* onePDMMx(){ return onePDMMx_.get();};
    inline TMap* onePDMMy(){ return onePDMMy_.get();};
    inline TMap* onePDMMz(){ return onePDMMz_.get();};

//  inline bool     isScattered(){ return isScattered_;};
    inline CQMemManager * memManager()     { return this->memManager_;     };

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


    void rotateDensities(const std::array<double,3>&,double);

/*
    template<typename OpIn, typename OpOut>
    static void ReImSeparate(OpIn&, 
        std::vector<std::reference_wrapper<OpOut>>&);

    template<typename OpIn, typename OpOut>
    static void ReImCombine(OpOut&, 
        std::vector<std::reference_wrapper<OpIn>>&);

    void sepReImOnePDM();
    void comReImOnePDM();
*/


    // MPI Routines
//  void mpiBCastDensity();
  };

/*
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
*/

  #include <quantum/quantum_scattergather.h>
  #include <quantum/quantum_misc.h>


/*
  template<>
  template<typename OpIn, typename OpOut>
  void Quantum<dcomplex>::ReImSeparate(OpIn &op, 
      std::vector<std::reference_wrapper<OpOut>> &sep){

    sep[0].get() = op.real();
    sep[1].get() = op.imag();
  };

  template<>
  template<typename OpIn, typename OpOut>
  void Quantum<double>::ReImSeparate(OpIn &op, 
      std::vector<std::reference_wrapper<OpOut>> &sep){;};


  template<>
  template<typename OpIn, typename OpOut>
  void Quantum<dcomplex>::ReImCombine(OpOut &op, 
      std::vector<std::reference_wrapper<OpIn>> &sep) {
  
      op.real() = sep[0].get();
      op.imag() = sep[1].get();
  };
*/

  /*
  template<>
  template<typename OpIn, typename OpOut>
  void Quantum<double>::ReImCombine(OpOut &op, 
      std::vector<std::reference_wrapper<OpIn>> &sep) {; };


  template<>
  inline void Quantum<double>::sepReImOnePDM() {};
  template<>
  inline void Quantum<double>::comReImOnePDM() {};

  template<>
  inline void Quantum<dcomplex>::sepReImOnePDM() {
    this->scatterDensity();
    std::vector<std::reference_wrapper<RealMatrix>> sep;

    // Separate Scalar OnePDM
    sep.emplace_back(*this->ReOnePDMScalar_);
    sep.emplace_back(*this->ImOnePDMScalar_);
    if(this->nTCS_ == 1 && this->isClosedShell)
      ReImSeparate(*this->onePDMA_,sep);
    else
      ReImSeparate(*this->onePDMScalar_,sep);

    sep.clear();

    if(this->nTCS_ == 2 || !this->isClosedShell) {
      sep.emplace_back(*this->ReOnePDMMz_);
      sep.emplace_back(*this->ImOnePDMMz_);
      ReImSeparate(*this->onePDMMz_,sep);
      sep.clear();
    }

    if(this->nTCS_ == 2) {
      sep.emplace_back(*this->ReOnePDMMx_);
      sep.emplace_back(*this->ImOnePDMMx_);
      ReImSeparate(*this->onePDMMx_,sep);
      sep.clear();

      sep.emplace_back(*this->ReOnePDMMy_);
      sep.emplace_back(*this->ImOnePDMMy_);
      ReImSeparate(*this->onePDMMy_,sep);
      sep.clear();
    };
  };

  template<>
  inline void Quantum<dcomplex>::comReImOnePDM() {
    this->scatterDensity();
    std::vector<std::reference_wrapper<RealMatrix>> sep;

    // Separate Scalar OnePDM
    sep.emplace_back(*this->ReOnePDMScalar_);
    sep.emplace_back(*this->ImOnePDMScalar_);
    if(this->nTCS_ == 1 && this->isClosedShell)
      ReImCombine(*this->onePDMA_,sep);
    else
      ReImCombine(*this->onePDMScalar_,sep);

    sep.clear();

    if(this->nTCS_ == 2 || !this->isClosedShell) {
      sep.emplace_back(*this->ReOnePDMMz_);
      sep.emplace_back(*this->ImOnePDMMz_);
      ReImCombine(*this->onePDMMz_,sep);
      sep.clear();
    }

    if(this->nTCS_ == 2) {
      sep.emplace_back(*this->ReOnePDMMx_);
      sep.emplace_back(*this->ImOnePDMMx_);
      ReImCombine(*this->onePDMMx_,sep);
      sep.clear();

      sep.emplace_back(*this->ReOnePDMMy_);
      sep.emplace_back(*this->ImOnePDMMy_);
      ReImCombine(*this->onePDMMy_,sep);
      sep.clear();
    };
  };
  */


}; // namespace ChronusQ



#endif
