
template<typename F>
class DIIS {
  // Useful Typedefs
  typedef Eigen::Matrix<F,Dynamic,Dynamic,ColMajor> FMat;
  typedef Eigen::Map<FMat> FMap;

  std::size_t               nExtrap_;
  std::size_t               eDim_;
  std::vector<F>            coeffs_;
  std::vector<H5::DataSet*> eFiles_;
  const std::function<void(H5::DataSet*,std::size_t,F*)> &read_;
  const std::function<F*(std::size_t)> &allocF_;
  const std::function<void(F*,std::size_t)> &freeF_;

  bool doNeg_;

public:

  DIIS() = delete;
  DIIS(std::size_t nExtrap, std::size_t eDim,
    const std::function<void(H5::DataSet*,std::size_t,F*)> &read,
    const std::function<F*(std::size_t)> &allocF,
    const std::function<void(F*,std::size_t)> &freeF,
    std::vector<H5::DataSet*> &eFiles, bool doNeg = true):
    nExtrap_(nExtrap), eDim_(eDim), read_(read), 
    eFiles_(eFiles), allocF_(allocF),freeF_(freeF),
    doNeg_(doNeg){

    coeffs_.resize(nExtrap_+1);
  }

  bool extrapolate();


  F* coeffs() { return coeffs_.data(); }
};

template<typename F>
bool DIIS<F>::extrapolate(){

  int N(nExtrap_ + 1);
  int NRHS(1); 
  int INFO(1);
  int LWORK(5*N);
  char NORM = 'O';
  double RCOND(0);
  bool InvFail(false);

  std::vector<F> BMem(N*N);
  FMap B(BMem.data(),N,N);
  std::vector<int> iPiv(N), iWork(N);

  F* SCRATCH1 = allocF_(eDim_);
  F* SCRATCH2 = allocF_(eDim_);

  FMap S1(SCRATCH1,eDim_,1);
  FMap S2(SCRATCH2,eDim_,1);

  B.setZero();
  for(auto iFile : eFiles_) {
    for(auto j = 0ul; j < nExtrap_; j++){
      read_(iFile,j,SCRATCH1);
    for(auto k = 0ul; k <= j; k++){
      read_(iFile,k,SCRATCH2);
      B(j,k) += S1.frobInner(S2);
    }
    }
  }

  B = B.template selfadjointView<Lower>();

  double fact = -1.0;
  if(!doNeg_) fact = 1.0;

  for(auto l = 0; l < nExtrap_; l++){
     B(nExtrap_,l) = fact;
     B(l,nExtrap_) = fact;
  }
  B(nExtrap_,nExtrap_) = 0.0;

  std::fill(coeffs_.begin(),coeffs_.end(),0.0);
  coeffs_[nExtrap_] = fact;

  prettyPrint(cout,B,"B");

  double ANORM = B.template lpNorm<1>();
  std::vector<F> WORK(LWORK);

  if(std::is_same<F,dcomplex>::value) {
    std::vector<double> RWORK(LWORK);
    // Linear Solve
    zgesv_(&N,&NRHS,reinterpret_cast<dcomplex*>(B.data()),&N,
        iPiv.data(),reinterpret_cast<dcomplex*>(coeffs_.data()),&N,&INFO);

    InvFail = (INFO != 0);

    // Obtain condition number from LU given from ZGESV
    zgecon_(&NORM,&N,reinterpret_cast<dcomplex*>(B.data()),&N,
        &ANORM,&RCOND,reinterpret_cast<dcomplex*>(WORK.data()),RWORK.data(),
        &INFO);

  } else if(std::is_same<F,double>::value) {
    // Linear Solve
    dgesv_(&N,&NRHS,reinterpret_cast<double*>(B.data()),&N,
        iPiv.data(),reinterpret_cast<double*>(coeffs_.data()),&N,&INFO);

    InvFail = (INFO != 0);

    // Obtain condition number from LU given from DGESV
    dgecon_(&NORM,&N,reinterpret_cast<double*>(B.data()),&N,
        &ANORM,&RCOND,reinterpret_cast<double*>(WORK.data()),iWork.data(),
        &INFO);
  }

  freeF_(SCRATCH1,eDim_);
  freeF_(SCRATCH2,eDim_);

  return !InvFail;
};

