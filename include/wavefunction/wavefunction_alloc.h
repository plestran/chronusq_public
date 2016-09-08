template <typename T>
void WaveFunction<T>::alloc() {
  this->checkMeta();
  Quantum<T>::alloc(this->nBasis_);

  auto NB = this->nBasis_;
  auto NBT = this->nTCS_ * NB;
  auto NBSq = NB*NB;
  auto NBTSq = NBT*NBT;

  // MO Coeffs
  this->moA_ = std::unique_ptr<TMap>(
    new TMap(this->memManager_->template malloc<T>(NBTSq),NBT,NBT)); 

  this->moA_->setZero();
  if(this->nTCS_ == 1 and !this->isClosedShell) {
    this->moB_ = std::unique_ptr<TMap>(
      new TMap(this->memManager_->template malloc<T>(NBTSq),NBT,NBT)); 

    this->moB_->setZero();
  }

  // MO Eigenenergies
  this->epsA_ = std::unique_ptr<RealMap>(
    new RealMap(this->memManager_->template malloc<double>(NBTSq),NBT,NBT)); 
    this->epsA_->setZero();
  if(this->nTCS_ == 1 and !this->isClosedShell){
    this->epsB_ = std::unique_ptr<RealMap>(
      new RealMap(this->memManager_->template malloc<double>(NBTSq),NBT,NBT)); 
    this->epsB_->setZero();
  }
}
