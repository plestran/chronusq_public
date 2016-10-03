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
  bool isAllocated_;

  CQMemManager * memManager_;


  // Pointers to TMatrix quantities that will store the density of
  // the quantum system in a finite basis

  // The Density in the basis of the identity and the Pauli matricies 
  std::unique_ptr<TMap> onePDMScalar_;
  std::unique_ptr<TMap> onePDMMz_;
  std::unique_ptr<TMap> onePDMMy_;
  std::unique_ptr<TMap> onePDMMx_;

  std::vector<TMap*> onePDM_;
  


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
    

  template<typename Scalar, DENSITY_TYPE DenTyp, typename Op>
  Scalar OperatorSpinCombine(const Op& op) {
    double zero = 0.0;
    bool isReal = typeid(T).hash_code() == typeid(double).hash_code();

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
  }

  public:

  bool isClosedShell;

  Quantum(){
    this->isClosedShell = false;
    this->nTCS_ = 1;
    this->maxMultipole_ = 3;
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
    memManager_(other.memManager_) {

    this->isAllocated_ = false;
    auto NB = other.onePDMScalar_->rows(); 
    this->alloc(NB);
   
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


    this->onePDMScalar_ = std::unique_ptr<TMap>(
        new TMap(this->memManager_->template malloc<T>(NSq),N,N));

    this->onePDMScalar_->setZero();
    this->onePDM_.emplace_back(this->onePDMScalar_.get());


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


  inline void setMaxMultipole(int i){ this->maxMultipole_ = i;   };
  inline void setNTCS(int i){         this->nTCS_ = i;           };

  inline int   nTCS(){ return nTCS_;};      
  inline int maxMultipole(){ return maxMultipole_;};

  inline std::vector<TMap*>& onePDM(){ return this->onePDM_;};
  inline TMap* onePDMScalar(){ return onePDMScalar_.get();};
  inline TMap* onePDMMx(){ return onePDMMx_.get();};
  inline TMap* onePDMMy(){ return onePDMMy_.get();};
  inline TMap* onePDMMz(){ return onePDMMz_.get();};

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



  // MPI Routines
//void mpiBCastDensity();
};
}; // namespace ChronusQ
