template <typename T>
void QuasiNewton2<T>::iniScratchFiles(){
  (*this->out_) << "Initializing Files for QuasiNewton Calculation" << endl;

  cout << "MAXSS " << this->maxSubSpace_ << endl;
  auto N      = this->qnObj_->nSingleDim();
  std::vector<hsize_t> dims;
  dims.push_back(this->maxSubSpace_);
  dims.push_back(N);

  std::vector<hsize_t> dimsSmall;
  dimsSmall.push_back(2*this->maxSubSpace_);
  dimsSmall.push_back(2*this->maxSubSpace_);

  time_t t; time(&t); // Time flag for files to be written

  std::string TRName = "Trial Vectors (Right) "    + std::to_string(t);
  std::string SRName = "Sigma Vectors (Right) "    + std::to_string(t);
  std::string PRName = "Rho Vectors (Right) "      + std::to_string(t);
  std::string RRName = "Residual Vectors (Right) " + std::to_string(t);
  std::string URName = "Solution Vectors (Right) " + std::to_string(t);
  std::string ASName = "A Super Matrix "           + std::to_string(t);
  std::string SSName = "S Super Matrix "           + std::to_string(t);

  this->TRFile_     = this->genScrFile_(H5PredType<T>(),TRName,dims);
//this->SigmaRFile_ = this->genScrFile_(H5PredType<T>(),SRName,dims);
//this->ResRFile_   = this->genScrFile_(H5PredType<T>(),RRName,dims);
//this->URFile_     = this->genScrFile_(H5PredType<T>(),URName,dims);

//if(this->matrixType_ == HERMETIAN_GEP)
//  this->RhoRFile_  = this->genScrFile_(H5PredType<T>(),PRName,dims);



  if(this->qnObj_->needsLeft_ or this->qnObj_->specialAlgorithm_ == SYMMETRIZED_TRIAL){
    std::string TLName = "Trial Vectors (Left) "    + std::to_string(t);
    std::string SLName = "Sigma Vectors (Left) "    + std::to_string(t);
    std::string PLName = "Rho Vectors (Left) "      + std::to_string(t);
    std::string RLName = "Residual Vectors (Left) " + std::to_string(t);
    std::string ULName = "Solution Vectors (Left) " + std::to_string(t);
   
    this->TLFile_     = this->genScrFile_(H5PredType<T>(),TLName,dims);
  //this->SigmaLFile_ = this->genScrFile_(H5PredType<T>(),SLName,dims);
  //this->ResLFile_   = this->genScrFile_(H5PredType<T>(),RLName,dims);
  //this->ULFile_     = this->genScrFile_(H5PredType<T>(),ULName,dims);
   
  //if(this->matrixType_ == HERMETIAN_GEP)
  //  this->RhoLFile_  = this->genScrFile_(H5PredType<T>(),PLName,dims);
  }

  if(this->qnObj_->specialAlgorithm_ == SYMMETRIZED_TRIAL) {
    this->ASuperFile_ = this->genScrFile_(H5PredType<T>(),ASName,dimsSmall);
    this->SSuperFile_ = this->genScrFile_(H5PredType<T>(),SSName,dimsSmall);
  }
}; // QuasiNewton2<T>::iniScratchFiles

template <typename T>
void QuasiNewton2<T>::writeTrialVectors(const int NTrial){
  auto N      = this->qnObj_->nSingleDim();
  hsize_t offset[] = {0,0};
  hsize_t stride[] = {1,1};
  hsize_t block[]  = {1,1};

  hsize_t subDim[] = {NTrial,N};
  hsize_t count[]  = {NTrial,N};

  cout << N << " " << NTrial << endl;
  H5::DataSpace memSpace(2,subDim,NULL);
  H5::DataSpace subDataSpace = this->TRFile_->getSpace();
  subDataSpace.selectHyperslab(H5S_SELECT_SET,count,offset,stride,block);
  this->TRFile_->write(this->TRMem_,H5PredType<T>(),memSpace,subDataSpace);
  if(this->qnObj_->needsLeft_ or this->qnObj_->specialAlgorithm_ == SYMMETRIZED_TRIAL)
    this->TLFile_->write(this->TLMem_,H5PredType<T>(),memSpace,subDataSpace);

}; // QuasiNewton2<T>::writeTrialVectors

template <typename T>
void QuasiNewton2<T>::readTrialVectors(const int NTrial){
  auto N      = this->qnObj_->nSingleDim();
  hsize_t offset[] = {0,0};
  hsize_t stride[] = {1,1};
  hsize_t block[]  = {1,1};

  hsize_t subDim[] = {NTrial,N};
  hsize_t count[]  = {NTrial,N};

  H5::DataSpace memSpace(2,subDim,NULL);
  H5::DataSpace subDataSpace = this->TRFile_->getSpace();
  subDataSpace.selectHyperslab(H5S_SELECT_SET,count,offset,stride,block);
  this->TRFile_->read(this->TRMem_,H5PredType<T>(),memSpace,
    subDataSpace);
  if(this->qnObj_->needsLeft_ or this->qnObj_->specialAlgorithm_ == SYMMETRIZED_TRIAL)
    this->TLFile_->read(this->TLMem_,H5PredType<T>(),memSpace,
      subDataSpace);

}; // QuasiNewton2<T>::readTrialVectors

