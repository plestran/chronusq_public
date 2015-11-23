template<typename T>
void Response<T>::initMeta() {
  this->checkWorkers();
  this->initPSCFDims();
  
  // Check that a method has been set
  if(this->iMeth_ == NOMETHOD)
    CErr("Must set a method to initialize Response Object",this->fileio_->out);

  // Determine type of Response Function
  if(this->iMeth_ == CIS || this->iMeth_ == RPA || this->iMeth_ == STAB)
    this->iClass_ = RESPONSE_CLASS::FOPPA;
  else if(this->iMeth_ == PPRPA || this->iMeth_ == PPTDA)
    this->iClass_ = RESPONSE_CLASS::PPPA;

  // Check if TDA
  if(this->iMeth_ == CIS || this->iMeth_ == PPTDA)
    this->doTDA_ = true;


  if(this->iClass_ == RESPONSE_CLASS::FOPPA) 
    this->initMetaFOPPA();
  else if(this->iClass_ == RESPONSE_CLASS::PPPA)
    this->initMetaPPRPA();

} // initMeta

template<typename T>
void Response<T>::initPSCFDims(){
  this->nBasis_         = this->singleSlater_->nBasis();
  this->nTCS_           = this->singleSlater_->nTCS();
  this->Ref_            = this->singleSlater_->Ref();
  this->nOA_            = this->singleSlater_->nOccA();
  this->nOB_            = this->singleSlater_->nOccB();
  this->nVA_            = this->singleSlater_->nVirA();
  this->nVB_            = this->singleSlater_->nVirB();

  this->nOAVA_          = this->nOA_*this->nVA_;
  this->nOBVB_          = this->nOB_*this->nVB_;
  this->nOAVB_          = this->nOA_*this->nVB_;
  this->nOBVA_          = this->nOB_*this->nVA_;
  this->nVAVA_SLT_      = this->nVA_*(this->nVA_-1)/2;
  this->nVBVB_SLT_      = this->nVB_*(this->nVB_-1)/2;
  this->nVAVA_LT_       = this->nVA_*(this->nVA_+1)/2;
  this->nVAVA_          = this->nVA_*this->nVA_;
  this->nVBVB_          = this->nVB_*this->nVB_;
  this->nOAOA_SLT_      = this->nOA_*(this->nOA_-1)/2;
  this->nOBOB_SLT_      = this->nOB_*(this->nOB_-1)/2;
  this->nOAOA_LT_       = this->nOA_*(this->nOA_+1)/2;
  this->nOAOA_          = this->nOA_*this->nOA_;
  this->nOBOB_          = this->nOB_*this->nOB_;
  this->nVAVB_          = this->nVA_*this->nVB_;
  this->nOAOB_          = this->nOA_*this->nOB_;
  this->nO_             = this->nOA_ + this->nOB_;
  this->nV_             = this->nVA_ + this->nVB_;
  this->nOV_            = this->nO_  * this->nV_;
  this->nVV_SLT_        = this->nV_*(this->nV_-1)/2;
  this->nVV_LT_         = this->nV_*(this->nV_+1)/2;
  this->nVV_            = this->nV_*this->nV_;
  this->nOO_SLT_        = this->nO_*(this->nO_-1)/2;
  this->nOO_LT_         = this->nO_*(this->nO_+1)/2;
  this->nOO_            = this->nO_*this->nO_;
}; // initPSCFDims

template<typename T>
void Response<T>::initMetaFOPPA(){

  if(this->iPart_ == SPIN_SEPARATED) {
    this->iMatIter_.push_back(FULL);
    if(this->Ref_ == SingleSlater<T>::TCS) {
      if(this->doTDA_) this->nMatDim_.push_back(this->nOV_  );
      else            this->nMatDim_.push_back(2*this->nOV_);
    } else {
      if(this->doTDA_) this->nMatDim_.push_back(this->nOAVA_  + this->nOBVB_  );
      else            this->nMatDim_.push_back(2*this->nOAVA_ + 2*this->nOBVB_);
    }
  } else if(iPart_ == SPIN_ADAPTED) {
    CErr("Spin-Adapted First Order Polarization Propagator NYI",
         this->fileio_->out);
  }

}; // initMetaFOPPA

template<typename T>
void Response<T>::initMetaPPRPA(){

  if(this->iPart_ == SPIN_SEPARATED) {
    if(this->Ref_ == SingleSlater<T>::TCS) {
      if(this->doTDA_) {
        this->iMatIter_.push_back(FULL_A_PPTDA);
        this->iMatIter_.push_back(FULL_C_PPTDA);
        this->nMatDim_.push_back(this->nVV_SLT_);
        this->nMatDim_.push_back(this->nOO_SLT_);
      } else {
        this->iMatIter_.push_back(FULL);
        this->nMatDim_.push_back(this->nVV_SLT_ + this->nOO_SLT_);
      }
    } else {
      if(this->doTDA_) {
        if(this->doAllAlpha_) {
          this->iMatIter_.push_back(AAA_PPTDA);
          this->iMatIter_.push_back(CAA_PPTDA);
          this->nMatDim_.push_back(this->nVAVA_SLT_);
          this->nMatDim_.push_back(this->nOAOA_SLT_);
        }
        if(this->doMixedAB_) {
          this->iMatIter_.push_back(AAB_PPTDA);
          this->iMatIter_.push_back(CAB_PPTDA);
          this->nMatDim_.push_back(this->nVAVB_);
          this->nMatDim_.push_back(this->nOAOB_);
        }
        if(this->doAllBeta_ && !this->singleSlater_->isClosedShell) {
          this->iMatIter_.push_back(ABB_PPTDA);
          this->iMatIter_.push_back(CBB_PPTDA);
          this->nMatDim_.push_back(this->nVBVB_SLT_);
          this->nMatDim_.push_back(this->nOBOB_SLT_);
        }
      } else {
        if(this->doAllAlpha_) {
          this->iMatIter_.push_back(AA_PPRPA);
          this->nMatDim_.push_back(this->nVAVA_SLT_ + this->nOAOA_SLT_);
        }
        if(this->doMixedAB_) {
          this->iMatIter_.push_back(AB_PPRPA);
          this->nMatDim_.push_back(this->nVAVB_ + this->nOAOB_);
        }
        if(this->doAllBeta_ && !this->singleSlater_->isClosedShell) {
          this->iMatIter_.push_back(BB_PPRPA);
          this->nMatDim_.push_back(this->nVBVB_SLT_ + this->nOBOB_SLT_);
        }
      }
    }
  } else if(this->iPart_ == SPIN_ADAPTED){
    CErr("Spin-Adapted Particle-Particle Propagator NYI",this->fileio_->out);
  }
};
