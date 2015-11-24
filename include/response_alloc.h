template<typename T>
void Response<T>::alloc(){

  // If we're doing full, allocate space for all states
  if(this->doFull_){
    for(auto iDim = this->nMatDim_.begin();
        iDim != this->nMatDim_.end();
        iDim++) {
      this->transDen_.push_back(TMat(*iDim,*iDim));
      this->frequencies_.push_back(VectorXd(*iDim));
    }
  } else {
    for(auto iDim = this->nMatDim_.begin();
        iDim != this->nMatDim_.end();
        iDim++) {
      this->transDen_.push_back(TMat(*iDim,this->nSek_));
      this->frequencies_.push_back(VectorXd(this->nSek_));
    }
  }
};
