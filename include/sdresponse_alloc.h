template<typename T>
void SDResponse<T>::alloc(){
  this->checkValid();
  this->omega_       = std::unique_ptr<VectorXd>(new VectorXd(this->nSek_));
  this->transDen_    = std::unique_ptr<TCMMatrix>(new TCMMatrix(
                       this->nSingleDim_, this->nSek_));
  this->oscStrength_ = std::unique_ptr<RealMatrix>(new RealMatrix(
                       this->nSek_+1, this->nSek_+1));
  this->transDipole_ = std::unique_ptr<RealTensor3d>(new RealTensor3d(
                       this->nSek_+1, this->nSek_+1, 3));
}
