#ifndef INCLUDED_RESP_DECL
#define INCLUDED_RESP_DECL

#include <pscf.h>
#include <qn.h>

namespace ChronusQ {

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
public:
  ResponseMatrix() : QNCallable<T>(),
    pscf_(NULL)
    { ; };
  ResponseMatrix(PostSCF<T> *pscf): QNCallable<T>(),
    pscf_(pscf)
    { this->memManager_ = pscf_->memManager(); };

  virtual void formFull() = 0;
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

  std::vector<ResponseMatrix<T>*> iMat_;

public:
  
  Response() : 
    PostSCF<T>(),
    useIncoreInts_(false),
    doFull_       (false),
    debugIter_    (false),
    doTDA_        (false),
    iJob_         (QNProblemType::DIAGONALIZATION)
    { cout << "In Response Constructor" << endl; };

  Response(Response &other) :
    PostSCF<T>(dynamic_cast<PostSCF<T>&>(other))
    { ; };

  // Quantum compliant
  virtual void formDensity() = 0;
  virtual void computeSSq() = 0;

/*
  // QNCallable compliant
  virtual void linearTrans(TMap &,TMap &,TMap &,TMap &,TMap &,TMap &) = 0;
  virtual void formGuess() = 0;
  virtual void formDiag()  = 0;
*/

  inline void communicate(WaveFunction<T> &wfn, CQMemManager &memManager) {
    PostSCF<T>::communicate(wfn,memManager);
    cout << "In Response Communicate" << endl;
  };

  inline void initMeta(){
    PostSCF<T>::initMeta();
    cout << "In Response initMeta" << endl;
  }
};
};

#endif

