/***********************
 * Form Density Matrix *
 ***********************/
template<typename T>
void SingleSlater<T>::formDensity(){
  if(!this->haveMO)
    CErr("No MO coefficients available to form one-particle density matrix!",
         this->fileio_->out);

  this->densityA_->setZero();
  *this->densityA_ = this->moA_->block(0,0,this->nBasis_,this->nOccA_)*
                   this->moA_->block(0,0,this->nBasis_,this->nOccA_).adjoint();
  if(this->RHF_) *this->densityA_ *= math.two;
  else {
    *this->densityB_ = this->moB_->block(0,0,this->nBasis_,this->nOccB_)*
                   this->moB_->block(0,0,this->nBasis_,this->nOccB_).adjoint();
  }
  if(this->controls_->printLevel>=2) {
    prettyPrint(this->fileio_->out,(*this->densityA_),"Alpha Density");
    if(!this->RHF_) prettyPrint(this->fileio_->out,(*this->densityB_),"Beta Density");
  };
  this->haveDensity = true;
}
/************************
 * Compute Total Energy *
 ************************/
template<typename T>
void SingleSlater<T>::computeEnergy(){
  this->energyOneE = (*this->aointegrals_->oneE_).frobInner(this->densityA_->conjugate());
#ifndef USE_LIBINT
  this->energyTwoE = ((*this->coulombA_)-(*this->exchangeA_)).frobInner(this->densityA->conjugate());
#else
  this->energyTwoE = (*this->PTA_).frobInner(this->densityA_->conjugate());
#endif
  this->totalEnergy= this->energyOneE + this->energyTwoE + this->energyNuclei;
  this->printEnergy();
};

