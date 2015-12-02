template<typename T>
void SDResponse<T>::alloc(){
  this->checkValid();
  this->omega_       = std::unique_ptr<VectorXd>(new VectorXd(this->nSek_));
  this->transDen_    = std::unique_ptr<TMatrix>(new TMatrix(
                       this->nSingleDim_, this->nSek_));
  this->oscStrength_ = std::unique_ptr<RealMatrix>(new RealMatrix(
                       this->nSek_+1, this->nSek_+1));
  this->transDipole_ = std::unique_ptr<RealTensor3d>(new RealTensor3d(
                       this->nSek_+1, this->nSek_+1, 3));
}

template<typename T>
void SDResponse<T>::iniSDResponse( Molecule * molecule, BasisSet * basisSet, 
  MOIntegrals<T> * mointegrals, FileIO * fileio, Controls * controls, 
  SingleSlater<T> * singleSlater) {

  this->communicate(*molecule,*basisSet,*singleSlater,*mointegrals,*fileio,
                    *controls);
  this->initMeta();
  this->setNSek(this->controls_->SDNSek);
  this->setMeth(this->controls_->SDMethod);
  this->initMeth();
  this->alloc();

}
