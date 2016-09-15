template <typename T>
void QuasiNewton2<T>::readGuess(){
  (*this->out_) << "Reading the Guess in QuasiNetwon" << endl;
  auto N    = this->qnObj_->nSingleDim();
  auto NVec = this->qnObj_->nGuess();

  H5::DataSpace dataspace = this->qnObj_->guessFile()->getSpace();
  this->qnObj_->guessFile()->read(this->TRMem_,H5PredType<T>(),
    dataspace,dataspace);

  if(this->qnObj_->needsLeft())
    this->qnObj_->guessFile()->read(this->TLMem_,H5PredType<T>(),
      dataspace,dataspace);
};

