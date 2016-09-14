#ifndef INCLUDED_POLARP
#define INCLUDED_POLARP

#include <response/response_decl.h>

namespace ChronusQ {
template <typename T>
class FOPPropagator : public ResponseMatrix<T> {
  // Useful Eigen Typedefs
  typedef Eigen::Matrix<T,Dynamic,Dynamic,ColMajor> TMat;
  typedef Eigen::Matrix<T,Dynamic,1> TVec;
  typedef Eigen::Map<TMat> TMap;
  typedef Eigen::Map<TVec> TVecMap;

public:

  FOPPropagator() : ResponseMatrix<T>(){ };
  FOPPropagator(PostSCF<T> *pscf): ResponseMatrix<T>(pscf){ };

  // QNCallable compliant
  void linearTrans(TMap &,TMap &,TMap &,TMap &,TMap &,TMap &){ };
  void formGuess(){ };
  void formDiag() { };

  void formFull();
};

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
    FOPPropagator<T> mat(this);
    this->iMat_.template push_back(dynamic_cast<ResponseMatrix<T>*>(&mat));
    mat.formFull();
  }
};


template <typename T>
void FOPPropagator<T>::formFull() {
  MOIntegrals<T> moints;
  moints.communicate(*this->pscf_->reference(),*this->pscf_->memManager());
  moints.initMeta();
  moints.formVOVO();
  moints.formVVOO();
  
  this->nSingleDim_ = this->pscf_->nOAVA() + this->pscf_->nOBVB();
  this->fullMatMem_ = 
    this->memManager_->template malloc<T>(this->nSingleDim_*this->nSingleDim_);
  std::fill_n(this->fullMatMem_,this->nSingleDim_*this->nSingleDim_,0.0);


  for(auto J = 0; J < this->pscf_->nOA(); J++)
  for(auto B = 0; B < this->pscf_->nVA(); B++)
  for(auto I = 0; I < this->pscf_->nOA(); I++) 
  for(auto A = 0; A < this->pscf_->nVA(); A++) {
    auto AI = A + I*this->pscf_->nVA();
    auto BJ = B + J*this->pscf_->nVA();
    auto AIBJ = AI + BJ*this->nSingleDim_;
    auto AIBJAAAA = 
      A + I*this->pscf_->nVA() + B*this->pscf_->nOAVA() + 
      J*this->pscf_->nVA()*this->pscf_->nOAVA();
    auto ABIJAAAA = 
      A + B*this->pscf_->nVA() + I*this->pscf_->nVAVA() + 
      J*this->pscf_->nOA()*this->pscf_->nVAVA();

    this->fullMatMem_[AIBJ] = moints.VOVOAAAA()[AIBJAAAA] - 
                        moints.VVOOAAAA()[ABIJAAAA];

    if(AI == BJ)
      this->fullMatMem_[AIBJ] += 
        (*this->pscf_->reference()->epsA())(A + this->pscf_->nOA()) -
        (*this->pscf_->reference()->epsA())(I);
  }

  for(auto J = 0; J < this->pscf_->nOB(); J++)
  for(auto B = 0; B < this->pscf_->nVB(); B++)
  for(auto I = 0; I < this->pscf_->nOA(); I++) 
  for(auto A = 0; A < this->pscf_->nVA(); A++) {
    auto AI = A + I*this->pscf_->nVA();
    auto BJ = B + J*this->pscf_->nVB() + this->pscf_->nOAVA();
    auto AIBJ = AI + BJ*this->nSingleDim_;
    auto AIBJAABB = 
      A + I*this->pscf_->nVA() + B*this->pscf_->nOAVA() + 
      J*this->pscf_->nVB()*this->pscf_->nOAVA();

    this->fullMatMem_[AIBJ] = moints.VOVOAABB()[AIBJAABB]; 
  }

  for(auto J = 0; J < this->pscf_->nOB(); J++)
  for(auto B = 0; B < this->pscf_->nVB(); B++)
  for(auto I = 0; I < this->pscf_->nOB(); I++) 
  for(auto A = 0; A < this->pscf_->nVB(); A++) {
    auto AI = A + I*this->pscf_->nVB() + this->pscf_->nOAVA();
    auto BJ = B + J*this->pscf_->nVB() + this->pscf_->nOAVA();
    auto AIBJ = AI + BJ*this->nSingleDim_;
    auto AIBJBBBB = 
      A + I*this->pscf_->nVB() + B*this->pscf_->nOBVB() + 
      J*this->pscf_->nVB()*this->pscf_->nOBVB();
    auto ABIJBBBB = 
      A + B*this->pscf_->nVB() + I*this->pscf_->nVBVB() + 
      J*this->pscf_->nOB()*this->pscf_->nVBVB();

    this->fullMatMem_[AIBJ] = moints.VOVOBBBB()[AIBJBBBB] - 
                              moints.VVOOBBBB()[ABIJBBBB];

    if(AI == BJ)
      this->fullMatMem_[AIBJ] += 
        (*this->pscf_->reference()->epsA())(A + this->pscf_->nOB()) -
        (*this->pscf_->reference()->epsA())(I);
  }

/*
  this->nSingleDim_ = this->pscf_->nOAVA();
  this->fullMatMem_ = 
    this->memManager_->template malloc<T>(this->nSingleDim_*this->nSingleDim_);

  for(auto J = 0; J < this->pscf_->nOA(); J++)
  for(auto B = 0; B < this->pscf_->nVA(); B++)
  for(auto I = 0; I < this->pscf_->nOA(); I++) 
  for(auto A = 0; A < this->pscf_->nVA(); A++) {
    auto AI = A + I*this->pscf_->nVA();
    auto BJ = B + J*this->pscf_->nVA();
    auto AIBJ = AI + BJ*this->nSingleDim_;
    auto AIBJAAAA = 
      A + I*this->pscf_->nVA() + B*this->pscf_->nOAVA() + 
      J*this->pscf_->nVA()*this->pscf_->nOAVA();
    auto ABIJAAAA = 
      A + B*this->pscf_->nVA() + I*this->pscf_->nVAVA() + 
      J*this->pscf_->nOA()*this->pscf_->nVAVA();

//  this->fullMatMem_[AIBJ] = 2*moints.VOVOAAAA()[AIBJAAAA] - 
//                      moints.VVOOAAAA()[ABIJAAAA];
    this->fullMatMem_[AIBJ] = - moints.VVOOAAAA()[ABIJAAAA];

    if(AI == BJ)
      this->fullMatMem_[AIBJ] += 
        (*this->pscf_->reference()->epsA())(A + this->pscf_->nOA()) -
        (*this->pscf_->reference()->epsA())(I);
  }
*/

  TMap Full(this->fullMatMem_,this->nSingleDim_,this->nSingleDim_);
  Full = Full.template selfadjointView<Upper>();
  prettyPrint(cout,Full,"Full");
  prettyPrint(cout,Full - Full.adjoint(),"Full");


  Eigen::SelfAdjointEigenSolver<TMat> es;
  es.compute(Full);
  VectorXd Eig = es.eigenvalues();
  prettyPrintSmart(cout,Eig,"Eig");
  prettyPrintSmart(cout,Eig*phys.eVPerHartree,"Eig");
}

};
#endif
