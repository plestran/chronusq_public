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

  bool doStab_;
public:

  FOPPropagator(PostSCF<T> *pscf, ResponseSettings sett, bool doStab = false):
    ResponseMatrix<T>(pscf,sett), doStab_(doStab){ };
  FOPPropagator() : ResponseMatrix<T>(){ };

  // QNCallable compliant
  void linearTrans(TMap &,TMap &,TMap &,TMap &,TMap &,TMap &);
  void formGuess();
  void formDiag();

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
  bool doStab_;
  
public:
  FOPPA(QNProblemType typ, RESPONSE_MATRIX_PARTITION part, bool doTDA, 
    bool doStab = false) : 
    Response<T>(typ,part,doTDA), doStab_(doStab){ };
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
  if(!this->sett_.doTDA) {
    this->nSingleDim_ *= 2;
    if(!this->doStab_) {
      this->matrixType_       = HERMETIAN_GEP;
      this->specialAlgorithm_ = SYMMETRIZED_TRIAL;
    }
  }
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
    /// Builds Spin-Separated A Matrix (TOP LEFT for !TDA)
      
    // A(ai,bj) (AAAA) = 
    //   d(ai,bj)(E(a) (A) - E(i) (A)) + (ai | bj) (AAAA) - (ab | ij) (AAAA)
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
 
    // A(ai,bj) (AABB) = (ai | bj) (AABB) 
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

    // A(ai,bj) (BBAA) gets taken care of by symmetrization at end
 
    // A(ai,bj) (BBBB) = 
    //   d(ai,bj)(E(a) (B) - E(i) (B)) + (ai | bj) (BBBB) - (ab | ij) (BBBB)
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

    // Build Top Right B for RPA
    if( not this->sett_.doTDA) {
      // B(ai,bj) (AAAA) = (ai | bj) (AAAA) - (aj | bi) (AAAA)
      for(auto J = 0; J < this->pscf_->nOA(); J++)
      for(auto B = 0; B < this->pscf_->nVA(); B++)
      for(auto I = 0; I < this->pscf_->nOA(); I++) 
      for(auto A = 0; A < this->pscf_->nVA(); A++) {
        auto AI = A + I*this->pscf_->nVA();
        auto BJ = B + J*this->pscf_->nVA() + this->nSingleDim_/2;
        auto AIBJ = AI + BJ*this->nSingleDim_;
        auto AIBJAAAA = 
          A + I*this->pscf_->nVA() + B*this->pscf_->nOAVA() + 
          J*this->pscf_->nVA()*this->pscf_->nOAVA();
        auto AJBIAAAA = 
          A + J*this->pscf_->nVA() + B*this->pscf_->nOAVA() + 
          I*this->pscf_->nVA()*this->pscf_->nOAVA();
 
        this->fullMatrix_[AIBJ] = moints.VOVOAAAA()[AIBJAAAA] - 
                            moints.VOVOAAAA()[AJBIAAAA];
      }

      // B(ai,bj) (AABB) = (ai | bj) (AABB) 
      for(auto J = 0; J < this->pscf_->nOB(); J++)
      for(auto B = 0; B < this->pscf_->nVB(); B++)
      for(auto I = 0; I < this->pscf_->nOA(); I++) 
      for(auto A = 0; A < this->pscf_->nVA(); A++) {
        auto AI = A + I*this->pscf_->nVA();
        auto BJ = B + J*this->pscf_->nVA() + this->pscf_->nOAVA() +
          this->nSingleDim_/2;
        auto AIBJ = AI + BJ*this->nSingleDim_;
        auto AIBJAABB = 
          A + I*this->pscf_->nVA() + B*this->pscf_->nOAVA() + 
          J*this->pscf_->nVB()*this->pscf_->nOAVA();
 
        this->fullMatrix_[AIBJ] = moints.VOVOAABB()[AIBJAABB];
      }

      // B(ai,bj) (BBAA) = B(bj,ai) (AABB)
      for(auto AI = this->pscf_->nOAVA(); AI < this->nSingleDim_/2; AI++)
      for(auto BJ = 0            ; BJ < this->pscf_->nOAVA()      ; BJ++) {
        auto AIBJ = AI + BJ*this->nSingleDim_ 
          + this->nSingleDim_*this->nSingleDim_/2;
        auto BJAI = BJ + AI*this->nSingleDim_
          + this->nSingleDim_*this->nSingleDim_/2;

        this->fullMatrix_[AIBJ] = this->fullMatrix_[BJAI];
      }

      // B(ai,bj) (BBBB) = (ai | bj) (BBBB) - (aj | bi) (BBBB)
      for(auto J = 0; J < this->pscf_->nOA(); J++)
      for(auto B = 0; B < this->pscf_->nVA(); B++)
      for(auto I = 0; I < this->pscf_->nOA(); I++) 
      for(auto A = 0; A < this->pscf_->nVA(); A++) {
        auto AI = A + I*this->pscf_->nVA() + this->pscf_->nOAVA();
        auto BJ = B + J*this->pscf_->nVA() + this->pscf_->nOAVA() +
          this->nSingleDim_/2;
        auto AIBJ = AI + BJ*this->nSingleDim_;
        auto AIBJBBBB = 
          A + I*this->pscf_->nVB() + B*this->pscf_->nOBVB() + 
          J*this->pscf_->nVB()*this->pscf_->nOBVB();
        auto AJBIBBBB = 
          A + J*this->pscf_->nVB() + B*this->pscf_->nOBVB() + 
          I*this->pscf_->nVB()*this->pscf_->nOBVB();
 
        this->fullMatrix_[AIBJ] = moints.VOVOBBBB()[AIBJBBBB] - 
                            moints.VOVOBBBB()[AJBIBBBB];
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

  // THESE ARE GENERAL TO ALL SPIN PARTITIONS
  if(not this->sett_.doTDA){
    // Copies A**H to bottom right for !TDA
    for(auto AI = 0; AI < this->nSingleDim_/2; AI++)
    for(auto BJ = 0; BJ < this->nSingleDim_/2; BJ++){
      auto AIBJTL = AI + this->nSingleDim_ * BJ;
      auto AIBJBR = (AI + this->nSingleDim_/2) + 
        this->nSingleDim_ * (BJ + this->nSingleDim_/2);

      this->fullMatrix_[AIBJBR] = std::conj(this->fullMatrix_[AIBJTL]);
    }


    // Bottom Left B gets taken care of my symmetrization

  } // Not TDA

  TMap Full(this->fullMatrix_,this->nSingleDim_,this->nSingleDim_);
  Full = Full.template selfadjointView<Upper>();
//prettyPrint(cout,Full,"Full");
//prettyPrint(cout,Full - Full.adjoint(),"Full");

  if(not this->sett_.doTDA and not this->doStab_){
  //Full.block(this->nSingleDim_/2,0,this->nSingleDim_/2,this->nSingleDim_) 
  //  *= -1;
  //this->matrixType_ = NON_HERMETIAN;
  }

/*
  TVec Eig = phys.eVPerHartree*Full.eigenvalues().real();
  std::sort(Eig.data(),Eig.data()+Eig.size());
  prettyPrintSmart(cout,Eig,"E");
*/
//CErr();
}


template <typename T>
void FOPPA<T>::alloc() {
  cout << "In FOPPA alloc" << endl;
  this->respMats_.emplace_back(std::make_shared<FOPPropagator<T>>(this,
    ResponseSettings{this->part_,true,false,this->doTDA_},this->doStab_));

  for(auto MAT : this->respMats_) {
    this->iMat_.template 
      push_back(dynamic_cast<ResponseMatrix<T>*>(MAT.get()));
  }

  Response<T>::alloc();
};

template <typename T>
void FOPPropagator<T>::formDiag() {
  cout << "In formDiag FOPPA" << endl;
  if(this->sett_.part == FULL) {
    if(this->pscf_->nTCS() == 2){ 
      for(auto I = 0, AI = 0; I < this->pscf_->nO(); I++      )
      for(auto A = 0        ; A < this->pscf_->nV(); A++, AI++) 
        this->diag_[AI] = 
          (*this->pscf_->reference()->epsA())(A + this->pscf_->nO()) -
          (*this->pscf_->reference()->epsA())(I);
    } else {
      CErr("Full FOPPA NOT IMPLEMENTED FOR NON 2C METHODS",
         this->pscf_->fileio()->out);
    }
  } else {
    for(auto I = 0, AI = 0; I < this->pscf_->nOA(); I++      )
    for(auto A = 0        ; A < this->pscf_->nVA(); A++, AI++) 
      this->diag_[AI] = 
        (*this->pscf_->reference()->epsA())(A + this->pscf_->nOA()) -
        (*this->pscf_->reference()->epsA())(I);

    if(this->sett_.part == SPIN_SEPARATED) {
      double * EPSB = this->pscf_->isClosedShell ?
        this->pscf_->reference()->epsA()->data() :
        this->pscf_->reference()->epsB()->data();
     
      for(auto I = 0, AI = this->pscf_->nOAVA(); I < this->pscf_->nOB(); I++)
      for(auto A = 0                     ; A < this->pscf_->nVB(); A++, AI++) 
        this->diag_[AI] = EPSB[A + this->pscf_->nOB()] - EPSB[I];
    }
  }

  // This is common to all spin partitions
  if(not this->sett_.doTDA)
    std::copy_n(this->diag_,this->nSingleDim_/2,
      this->diag_ + this->nSingleDim_/2);
}


template <typename T>
void FOPPropagator<T>::formGuess() {

  std::vector<hsize_t> dims;
  dims.push_back(this->nGuess_);
  dims.push_back(this->nSingleDim_);

  this->guessFile_ = 
    this->pscf_->fileio()->createScratchPartition(H5PredType<T>(),
      "FOPPA Response Guess",dims);

  // Initialize an index vector with increasing ints
  std::vector<int> indx(this->nSingleDim_/2,0);
  std::iota(indx.begin(),indx.end(),0);


  // Sort the index vector based on the diagonal of the RM
  // (uses Lambda expression)
  std::sort(indx.begin(),indx.end(),
    [&](const int& a, const int& b){
      return this->diag_[a] < this->diag_[b];
    }
  );

  
  for(auto iGuess = 0; iGuess < this->nGuess_; iGuess++) {
    T one = T(1.0);

  
//  HDF5 Stores things RowMajor...
    hsize_t offset[] = {iGuess,indx[iGuess]};
    hsize_t count[]  = {1,1};
    hsize_t stride[] = {1,1};
    hsize_t block[]  = {1,1};
    hsize_t subDim[] = {1,1};

    H5::DataSpace memSpace(2,subDim,NULL);
    H5::DataSpace subDataSpace = this->guessFile_->getSpace();
    subDataSpace.selectHyperslab(H5S_SELECT_SET,count,offset,stride,block);
    this->guessFile_->write(&one,H5PredType<T>(),memSpace,subDataSpace);
    
  }

/*
  T* GUESS = this->memManager_->template malloc<T>(this->nGuess_*this->nSingleDim_);

  this->guessFile_->read(GUESS,H5PredType<T>(),this->guessFile_->getSpace(),
    this->guessFile_->getSpace());
  TMap G(GUESS,this->nSingleDim_,this->nGuess_);
  prettyPrintSmart(cout,G,"Guess");
*/

};
template <typename T>
void FOPPropagator<T>::linearTrans(TMap &TR,TMap &TL,TMap &SR,TMap &SL,
   TMap &RR,TMap &RL) {
  TMap Full(this->fullMatrix_,this->nSingleDim_,this->nSingleDim_);
  SR = Full * TR;
  prettyPrintSmart(cout,SR,"Sigma R");

  if(not this->sett_.doTDA and not this->doStab_) {
    SL = Full * TL;
    prettyPrintSmart(cout,SL,"Sigma L");
    RR = TL;
    RL = TR;
    RR.block(this->nSingleDim_/2,0,this->nSingleDim_/2,RR.cols()) *= -1;
    RL.block(this->nSingleDim_/2,0,this->nSingleDim_/2,RL.cols()) *= -1;
    prettyPrintSmart(cout,RR,"Rho R");
    prettyPrintSmart(cout,RL,"Rho L");
  }
}

};
#endif
