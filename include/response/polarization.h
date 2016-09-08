#ifndef INCLUDED_POLARP
#define INCLUDED_POLARP

#include <response/response_decl.h>

namespace ChronusQ {
template <typename T>
class FOPPA : public Response<T> {
  // Useful Eigen Typedefs
  typedef Eigen::Matrix<T,Dynamic,Dynamic,ColMajor> TMat;
  typedef Eigen::Matrix<T,Dynamic,1> TVec;
  typedef Eigen::Map<TMat> TMap;
  typedef Eigen::Map<TVec> TVecMap;

public:
  FOPPA() : Response<T>()
    { cout << "In FOPPA Constructor" << endl; };
  FOPPA(FOPPA &other) :
    Response<T>(dynamic_cast<Response<T>&>(other))
    { ; };

  // Quantum compliant
  void formDensity(){ };
  void computeSSq(){ };

/*
  // QNCallable compliant
  void linearTrans(TMap &,TMap &,TMap &,TMap &,TMap &,TMap &){ };
  void formGuess(){ };
  void formDiag() { };
*/

  inline void initMeta() {
    Response<T>::initMeta();
    cout << "In FOPPA initMeta" << endl;
  }
};

template <typename T>
class FOPPAResponseMatrix : public ResponseMatrix<T> {
  // Useful Eigen Typedefs
  typedef Eigen::Matrix<T,Dynamic,Dynamic,ColMajor> TMat;
  typedef Eigen::Matrix<T,Dynamic,1> TVec;
  typedef Eigen::Map<TMat> TMap;
  typedef Eigen::Map<TVec> TVecMap;

public:
  // QNCallable compliant
  void linearTrans(TMap &,TMap &,TMap &,TMap &,TMap &,TMap &){ };
  void formGuess(){ };
  void formDiag() { };
};

};
#endif
