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
  FOPPropagator(PostSCF<T> *pscf, ResponseSettings sett):
    ResponseMatrix<T>(pscf,sett){ };

  // QNCallable compliant
  void linearTrans(TMap &,TMap &,TMap &,TMap &,TMap &,TMap &){ };
  void formGuess(){ };
  void formDiag() { };

  // ResponseMatrix Compliant
  void formFull();
  void initMeta();
};

template <typename T>
class FOPPA : public Response<T> {
  // Useful Eigen Typedefs
  typedef Eigen::Matrix<T,Dynamic,Dynamic,ColMajor> TMat;
  typedef Eigen::Matrix<T,Dynamic,1> TVec;
  typedef Eigen::Map<TMat> TMap;
  typedef Eigen::Map<TVec> TVecMap;


  std::vector<std::shared_ptr<FOPPropagator<T>>> respMats_;
  
public:
  FOPPA(QNProblemType typ, RESPONSE_MATRIX_PARTITION part, bool doTDA) : 
    Response<T>(typ,part,doTDA){ };
  FOPPA() : Response<T>(){ };
  FOPPA(FOPPA &other) : Response<T>(dynamic_cast<Response<T>&>(other)){ };

  // Quantum compliant
  void formDensity(){ };
  void computeSSq(){ };

  inline void initMeta() {
    Response<T>::initMeta();
    cout << "In FOPPA initMeta" << endl;
  }

  // Response compliant
//void runResponse();

  void alloc();
};

template <typename T>
void FOPPropagator<T>::initMeta() {
  cout << "In FOPPropagator initMeta" << endl;
  if(this->sett_.part == FULL) {
    this->nSingleDim_ = this->pscf_->nOV();

  } else if(this->sett_.part == SPIN_SEPARATED) {

    this->nSingleDim_ = this->pscf_->nOAVA() + this->pscf_->nOBVB();

  } else if(this->sett_.part == SPIN_ADAPTED) {

    this->nSingleDim_ = this->pscf_->nOAVA();

  }
  if(!this->sett_.doTDA) this->nSingleDim_ *= 2;
}

template <typename T>
void FOPPropagator<T>::formFull() {
  MOIntegrals<T> moints;
  moints.communicate(*this->pscf_->reference(),*this->pscf_->memManager());
  moints.initMeta();
  moints.formVOVO();
  moints.formVVOO();
  
  this->fullMatrix_ = 
    this->memManager_->template malloc<T>(this->nSingleDim_*this->nSingleDim_);
  std::fill_n(this->fullMatrix_,this->nSingleDim_*this->nSingleDim_,0.0);


  if(this->sett_.part == SPIN_SEPARATED) {
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
 
      this->fullMatrix_[AIBJ] = moints.VOVOAAAA()[AIBJAAAA] - 
                          moints.VVOOAAAA()[ABIJAAAA];
 
      if(AI == BJ)
        this->fullMatrix_[AIBJ] += 
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
 
      this->fullMatrix_[AIBJ] = moints.VOVOAABB()[AIBJAABB]; 
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
 
      this->fullMatrix_[AIBJ] = moints.VOVOBBBB()[AIBJBBBB] - 
                                moints.VVOOBBBB()[ABIJBBBB];
 
      if(AI == BJ){
        if(this->pscf_->nTCS() == 1 and !this->pscf_->isClosedShell)
          this->fullMatrix_[AIBJ] += 
            (*this->pscf_->reference()->epsB())(A + this->pscf_->nOB()) -
            (*this->pscf_->reference()->epsB())(I);
        else
          this->fullMatrix_[AIBJ] += 
            (*this->pscf_->reference()->epsA())(A + this->pscf_->nOB()) -
            (*this->pscf_->reference()->epsA())(I);
      }
    }
  } else if(this->sett_.part == SPIN_ADAPTED) {

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

      if(this->sett_.doSinglets)
        this->fullMatrix_[AIBJ] = 2*moints.VOVOAAAA()[AIBJAAAA] - 
                            moints.VVOOAAAA()[ABIJAAAA];
      else if(this->sett_.doTriplets)
        this->fullMatrix_[AIBJ] = - moints.VVOOAAAA()[ABIJAAAA];

      if(AI == BJ)
        this->fullMatrix_[AIBJ] += 
          (*this->pscf_->reference()->epsA())(A + this->pscf_->nOA()) -
          (*this->pscf_->reference()->epsA())(I);
    }
  }

  TMap Full(this->fullMatrix_,this->nSingleDim_,this->nSingleDim_);
  Full = Full.template selfadjointView<Upper>();
//prettyPrint(cout,Full,"Full");
//prettyPrint(cout,Full - Full.adjoint(),"Full");


//Eigen::SelfAdjointEigenSolver<TMat> es;
//es.compute(Full);
//VectorXd Eig = es.eigenvalues();
//prettyPrintSmart(cout,Eig,"Eig");
//prettyPrintSmart(cout,Eig*phys.eVPerHartree,"Eig");
}


template <typename T>
void FOPPA<T>::alloc() {
  cout << "In FOPPA alloc" << endl;
  this->respMats_.emplace_back(std::make_shared<FOPPropagator<T>>(this,
    ResponseSettings{this->part_,true,false,this->doTDA_}));

  for(auto MAT : this->respMats_) {
    this->iMat_.template 
      push_back(dynamic_cast<ResponseMatrix<T>*>(MAT.get()));
  }

  Response<T>::alloc();
};

/*
template <typename T>
void FOPPA<T>::runResponse() {
  FOPPropagator<T> mat(this,
     ResponseSettings{this->part_,true,false,this->doTDA_});

  this->iMat_.template push_back(dynamic_cast<ResponseMatrix<T>*>(&mat));
  this->iMat_[0]->initMeta();
  this->iMat_[0]->formFull();

  std::function<H5::DataSet*(const H5::CompType&,std::string&,
    std::vector<hsize_t>&)> fileFactory = 
      std::bind(&FileIO::createScratchPartition,this->fileio_,
      std::placeholders::_1,std::placeholders::_2,std::placeholders::_3);

  QuasiNewton2<T> qn(this->iMat_[0],this->memManager_,fileFactory);
  qn.setAlgorithm(FULL_SOLVE);
  qn.run();
}
*/

};
#endif
