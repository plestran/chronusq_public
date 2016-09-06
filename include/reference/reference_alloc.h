template <typename T>
void Reference<T>::alloc() {
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
}
