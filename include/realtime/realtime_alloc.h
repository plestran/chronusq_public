template <typename T>
void RealTime<T>::alloc() {
  ssPropagator_ = std::unique_ptr<SingleSlater<dcomplex>>(
    new SingleSlater<dcomplex>(groundState_));
  ssPropagator_->isNotPrimary();

//prettyPrintSmart(cout,*groundState_->onePDMScalar(),"GSPS");
//prettyPrintSmart(cout,*groundState_->fockScalar(),"GSFS");
//prettyPrintSmart(cout,*ssPropagator_->onePDMScalar(),"PS");
//prettyPrintSmart(cout,*ssPropagator_->fockScalar(),"FS");

//CErr();
  auto NB  = ssPropagator_->basisset()->nBasis();
  auto NBT = NB * ssPropagator_->nTCS();

  NBSqScratch_  = memManager_->template malloc<dcomplex>(NB*NB);
  UTransScalar_ = memManager_->template malloc<dcomplex>(NB*NB);
  POScalarSav_  = std::unique_ptr<ComplexMap>(
    new ComplexMap(memManager_->template malloc<dcomplex>(NB*NB),NB,NB));

  POSav_.emplace_back(POScalarSav_.get());

  if(ssPropagator_->nTCS() == 2 or !ssPropagator_->isClosedShell){
    UTransMz_= memManager_->template malloc<dcomplex>(NB*NB);
    POMzSav_  = std::unique_ptr<ComplexMap>(
      new ComplexMap(memManager_->template malloc<dcomplex>(NB*NB),NB,NB));

    POSav_.emplace_back(POMzSav_.get());
  }
  if(ssPropagator_->nTCS() == 2){ 
    NBTSqScratch_  = memManager_->template malloc<dcomplex>(NBT*NBT);
    NBTSqScratch2_ = memManager_->template malloc<dcomplex>(NBT*NBT);
    UTransMx_= memManager_->template malloc<dcomplex>(NB*NB);
    UTransMy_= memManager_->template malloc<dcomplex>(NB*NB);
    POMySav_  = std::unique_ptr<ComplexMap>(
      new ComplexMap(memManager_->template malloc<dcomplex>(NB*NB),NB,NB));
    POMxSav_  = std::unique_ptr<ComplexMap>(
      new ComplexMap(memManager_->template malloc<dcomplex>(NB*NB),NB,NB));

    POSav_.emplace_back(POMySav_.get());
    POSav_.emplace_back(POMxSav_.get());
  } else {
    NBTSqScratch_ = NBSqScratch_;
  }
}

