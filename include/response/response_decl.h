#ifndef INCLUDED_RESP_DECL
#define INCLUDED_RESP_DECL

#include <pscf.h>
#include <qn.h>

namespace ChronusQ {

enum RESPONSE_MATRIX_PARTITION {
  FULL,
  SPIN_SEPARATED,
  SPIN_ADAPTED
};

struct ResponseSettings {
  RESPONSE_MATRIX_PARTITION part;
  bool doSinglets;
  bool doTriplets;
  bool doTDA;
};

inline ResponseSettings DefaultResponseSettings() {
  return ResponseSettings{FULL,true,false,true};
};

template <typename T>
class ResponseMatrix : public QNCallable<T> {
  // Useful Eigen Typedefs
  typedef Eigen::Matrix<T,Dynamic,Dynamic,ColMajor> TMat;
  typedef Eigen::Matrix<T,Dynamic,1> TVec;
  typedef Eigen::Map<TMat> TMap;
  typedef Eigen::Map<TVec> TVecMap;

protected:
  CQMemManager *memManager_;
  PostSCF<T>   *pscf_;


  T * fullMatMem_;

  ResponseSettings sett_;

public:

  ResponseMatrix(PostSCF<T> *pscf, ResponseSettings sett): QNCallable<T>(),
    pscf_(pscf), sett_(sett) { 
    this->memManager_ = pscf_->memManager(); 
    this->checkMeta();
  };

  ResponseMatrix(PostSCF<T> *pscf): 
    ResponseMatrix(pscf,DefaultResponseSettings()){ };

  ResponseMatrix() : ResponseMatrix(NULL,DefaultResponseSettings()){ };


  virtual void formFull() = 0;

  virtual void initMeta() = 0;
  inline void checkMeta() {
    if(this->pscf_->nTCS() == 2 and this->sett_.part != FULL){
      cout << "Switching Matrix partitioning to FULL for 2C method" << endl;
      this->sett_.part = FULL;
    } else if(this->pscf_->nTCS() == 1 and 
              !this->pscf_->isClosedShell and 
              this->sett_.part == SPIN_ADAPTED) {
      cout << "Switching Matrix partitioning to SPIN_SEPARATED for ";
      cout << "Unresticted method" << endl;
      this->sett_.part = SPIN_SEPARATED;
    }
  }
};

template <typename T>
class Response : public PostSCF<T> {
  // Useful Eigen Typedefs
  typedef Eigen::Matrix<T,Dynamic,Dynamic,ColMajor> TMat;
  typedef Eigen::Matrix<T,Dynamic,1> TVec;
  typedef Eigen::Map<TMat> TMap;
  typedef Eigen::Map<TVec> TVecMap;

protected:

  /** Job Control for Response calculation **/
  bool useIncoreInts_;   ///< Use Incore Integrals
  bool doFull_;          ///< Do full matrix problem (in-core)
  bool debugIter_;       ///< Diagonalize full incore matrix iteratively
  bool doTDA_;           ///< Invoke TDA

  QNProblemType iJob_;   ///< Response Job Type
  RESPONSE_MATRIX_PARTITION part_;

  std::vector<ResponseMatrix<T>*> iMat_;

public:
  
  Response(QNProblemType typ, RESPONSE_MATRIX_PARTITION part, bool doTDA): 
    PostSCF<T>(),
    iJob_(typ), doTDA_(doTDA), part_(part),
    useIncoreInts_(false),
    doFull_       (false),
    debugIter_    (false){ };

  Response() : Response(QNProblemType::DIAGONALIZATION,false){ };

  Response(Response &other) :
    PostSCF<T>(dynamic_cast<PostSCF<T>&>(other))
    { ; };

  // Quantum compliant
  virtual void formDensity() = 0;
  virtual void computeSSq() = 0;


  inline void communicate(WaveFunction<T> &wfn, CQMemManager &memManager) {
    PostSCF<T>::communicate(wfn,memManager);
    cout << "In Response Communicate" << endl;
  };

  inline void initMeta(){
    PostSCF<T>::initMeta();
    cout << "In Response initMeta" << endl;
  }

  virtual void runResponse() = 0;


  void doTDA()  { this->doTDA_  = true; };
  void doFull() { this->doFull_ = true; };
};
};

#endif

