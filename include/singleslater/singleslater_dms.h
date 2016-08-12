template<typename T>
void SingleSlater<T>::formADMPGrad(int IDMSIter){
  prettyPrint(cout,(*this->fockOrthoA_),"FA"); 
  prettyPrint(cout,(*this->onePDMOrthoA_),"DA"); 
  this->NBSqScratch2_->noalias() =
    (*this->fockOrthoA_) * (*this->onePDMOrthoA_);

  // G = 3 * (FP + PF)
  this->NBSqScratch_->noalias() = 
    3 * ((*this->NBSqScratch2_) + this->NBSqScratch2_->adjoint());

  // G -= 2 ( FPP + PFP + PPF )
  this->NBSqScratch_->noalias() -= 
    2 * (*this->NBSqScratch2_) * (*this->onePDMOrthoA_); 
  this->NBSqScratch_->noalias() -= 
    2 * (*this->onePDMOrthoA_) * (*this->NBSqScratch2_); 
  this->NBSqScratch_->noalias() -= 
    2 * (*this->onePDMOrthoA_) * this->NBSqScratch2_->adjoint(); 

  this->writeDIIS(this->ADMPGradScalar_,IDMSIter % this->nDIISExtrap_,
    this->NBSqScratch_->data());  
};

template<typename T>
void SingleSlater<T>::formDMSErr(int IDMSIter){
  this->formADMPGrad(IDMSIter);

  this->readDIIS(this->ADMPGradScalar_,IDMSIter % this->nDIISExtrap_,
    this->NBSqScratch_->data());  
  prettyPrint(cout,(*this->NBSqScratch_),"Grad");

  this->NBSqScratch2_->noalias() =
    (*this->onePDMOrthoA_) * (*this->fockOrthoA_)  ;

  for(auto mu = 0; mu < this->nBasis_; mu++)
  for(auto nu = 0; nu < this->nBasis_; nu++) {
    (*this->NBSqScratch_)(mu,nu) /=
      (3.0 + 2.0 * (*this->onePDMA_)(mu,mu)) * (*this->fockA_)(nu,nu) +
      (3.0 + 2.0 * (*this->onePDMA_)(nu,nu)) * (*this->fockA_)(mu,mu) -
      4.0 * ((*this->NBSqScratch2_)(mu,mu) + (*this->NBSqScratch2_)(nu,nu));
  }
 
  (*this->NBSqScratch_) *= -1;
  prettyPrint(cout,(*this->NBSqScratch_),"DMS E");
  this->writeDIIS(this->DMSErrScalar_,IDMSIter % this->nDIISExtrap_,
    this->NBSqScratch_->data());  

};

template<typename T>
void SingleSlater<T>::DMSExtrap(int NDIIS) {
  std::function<void(H5::DataSet*,std::size_t,T*)> F1 = 
    std::bind(&SingleSlater<T>::readDIIS,this,
      std::placeholders::_1,std::placeholders::_2,std::placeholders::_3);

  std::function<T*(std::size_t)> F2 = 
    std::bind(&CQMemManager::malloc<T>,this->memManager_,
      std::placeholders::_1);

  std::function<void(T*,std::size_t)> F3 = 
    std::bind(&CQMemManager::free<T>,this->memManager_,
      std::placeholders::_1,std::placeholders::_2);

  T* coef;

  auto Extrap = [&](TMap *res, H5::DataSet *basis) {
    for(auto j = 0; j < NDIIS; j++) {
      this->readDIIS(basis,j,this->NBSqScratch_->data());
      (*res) += coef[j] * (*this->NBSqScratch_); 
    }
  };

  DIIS<T> extrap(NDIIS,this->nBasis_*this->nBasis_,F1,F2,F3,
    this->DMSErr_,false);

  
  if(!extrap.extrapolate()) {
    cout << "INVFailed" << endl;
  CErr();
    return;
  }
  coef = extrap.coeffs();
  cout << "COEFFS" << endl;
  for(auto i = 0; i < NDIIS; i++) cout <<coef[i] <<endl;
  if(NDIIS == this->nDIISExtrap_-3) CErr();

  for(auto DEN : this->onePDM_) DEN->setZero();

  for(auto I = 0; I < this->onePDM_.size(); I++){
    Extrap(this->onePDM_[I],this->DDIIS_[I]);
    Extrap(this->onePDM_[I],this->DMSErr_[I]);
  }
  prettyPrint(cout,(*this->onePDMA_),"P After Exp");

  
  for(auto I = 0; I < this->onePDM_.size(); I++){
    (*this->onePDM_[I]) *= 
      static_cast<T>(2*this->molecule_->nTotalE())/this->onePDM_[0]->trace();
  }
  prettyPrint(cout,(*this->onePDMA_),"P After Scale");
  
  this->McWeeny(10);
  prettyPrint(cout,(*this->onePDMA_),"P After Pure");
  CErr();
}

template<typename T>
void SingleSlater<T>::McWeeny(int MAX){
  for(auto I = 0; I < MAX; I++){
    prettyPrint(cout,(*this->onePDMA_),"P in McWeeny " + std::to_string(I));
    this->NBSqScratch_->noalias() = (*this->onePDMA_) * (*this->onePDMA_);
    prettyPrint(cout,(*this->NBSqScratch_),"P2 in McWeeny " + std::to_string(I));
    (*this->NBSqScratch2_) = (*this->onePDMA_) - (*this->NBSqScratch_);
    if(this->NBSqScratch2_->norm() < 1e-10) break;

    (*this->NBSqScratch2_) = (*this->onePDMA_);
    this->onePDMA_->noalias() = 3 * (*this->NBSqScratch_);
    this->onePDMA_->noalias() -= 
      2 * (*this->NBSqScratch_) * (*this->NBSqScratch2_);
  }
};

template<typename T>
void SingleSlater<T>::initDMSFiles(){
  if(!this->isPrimary) this->doDMS = false;
  if(!this->doDMS) return;

  this->initADMPFiles();
  std::vector<hsize_t> dims;
  dims.push_back(this->nDIISExtrap_);
  dims.push_back(this->nBasis_);
  dims.push_back(this->nBasis_);

  this->DMSErrScalar_ = 
    this->fileio_->createScratchPartition(H5PredType<T>(),
      "DMS Error (Scalar) For  Extrapoloation",dims);

  this->DMSErr_.emplace_back(this->DMSErrScalar_);

  if(this->nTCS_ == 2 || !this->isClosedShell) {
    this->DMSErrMz_ = 
      this->fileio_->createScratchPartition(H5PredType<T>(),
        "DMS Error (Mz) For  Extrapoloation",dims);
    this->DMSErr_.emplace_back(this->DMSErrMz_);
  }

  if(this->nTCS_ == 2) {
    this->DMSErrMy_ = 
      this->fileio_->createScratchPartition(H5PredType<T>(),
        "DMS Error (My) For  Extrapoloation",dims);

    this->DMSErrMx_ = 
      this->fileio_->createScratchPartition(H5PredType<T>(),
        "DMS Error (Mx) For  Extrapoloation",dims);
    this->DMSErr_.emplace_back(this->DMSErrMy_);
    this->DMSErr_.emplace_back(this->DMSErrMz_);
  }
}

template<typename T>
void SingleSlater<T>::initADMPFiles(){
  std::vector<hsize_t> dims;
  dims.push_back(this->nDIISExtrap_);
  dims.push_back(this->nBasis_);
  dims.push_back(this->nBasis_);

  this->ADMPGradScalar_ = 
    this->fileio_->createScratchPartition(H5PredType<T>(),
      "ADMP Gradient (Scalar) For  Extrapoloation",dims);

  if(this->nTCS_ == 2 || !this->isClosedShell) {
    this->ADMPGradMz_ = 
      this->fileio_->createScratchPartition(H5PredType<T>(),
        "ADMP Gradient (Mz) For  Extrapoloation",dims);
  }

  if(this->nTCS_ == 2) {
    this->ADMPGradMy_ = 
      this->fileio_->createScratchPartition(H5PredType<T>(),
        "ADMP Gradient (My) For  Extrapoloation",dims);

    this->ADMPGradMx_ = 
      this->fileio_->createScratchPartition(H5PredType<T>(),
        "ADMP Gradient (Mx) For  Extrapoloation",dims);
  }
}
