#include <response.h>

namespace ChronusQ {

/**
 *  Places the vectorized contents of T into the virtual-occupied block
 *  of TMOA and TMOB
 */
template<>
void Response<double>::placeVirOcc(RealVecMap &T, RealMatrix &TMOA,
  RealMatrix &TMOB) {

  bool doBeta = this->iPart_ != SPIN_ADAPTED && 
                this->Ref_ != SingleSlater<double>::TCS &&
                this->iClass_ != PPPA;

  Eigen::Map<RealMatrix>  TExpandedA(T.data(),this->nVA_,this->nOA_);
  Eigen::Map<RealMatrix>  TExpandedB(T.data(),0,0);
  if(doBeta)
    new (&TExpandedB) Eigen::Map<RealMatrix>(
      T.data()+this->nOAVA_,this->nVB_,this->nOB_
    );
  
  auto iOffA = this->nOA_;
  auto iOffB = this->nOB_;

  TMOA.block(iOffA,0,this->nVA_,this->nOA_) = TExpandedA;
  if(doBeta)
    TMOB.block(iOffB,0,this->nVA_,this->nOA_) = TExpandedB;

}; // placeVirOcc

template<>
void Response<double>::placeOccVir(RealVecMap &T, RealMatrix &TMOA,
  RealMatrix &TMOB) {

  bool doBeta = this->iPart_ != SPIN_ADAPTED && 
                this->Ref_ != SingleSlater<double>::TCS &&
                this->iClass_ != PPPA;

  Eigen::Map<RealMatrix>  TExpandedA(T.data(),this->nVA_,this->nOA_);
  Eigen::Map<RealMatrix>  TExpandedB(T.data(),0,0);
  if(doBeta)
    new (&TExpandedB) Eigen::Map<RealMatrix>(
      T.data()+this->nOAVA_,this->nVB_,this->nOB_
    );
  
  auto iOffA = this->nOA_;
  auto iOffB = this->nOB_;

  TMOA.block(0,iOffA,this->nOA_,this->nVA_) = TExpandedA.adjoint();
  if(doBeta)
    TMOB.block(0,iOffB,this->nOA_,this->nVA_) = TExpandedB.adjoint();

}; // placeVirOcc

template<>
void Response<double>::retrieveVirOcc(RealVecMap &T, RealMatrix &TMOA,
  RealMatrix &TMOB) {

  bool doBeta = this->iPart_ != SPIN_ADAPTED && 
                this->Ref_ != SingleSlater<double>::TCS &&
                this->iClass_ != PPPA;

  auto iOffA = this->nOA_;
  auto iOffB = this->nOB_;

  RealMatrix TExpandedA, TExpandedB;
  TExpandedA = TMOA.block(iOffA,0,this->nVA_,this->nOA_);
  if(doBeta) TExpandedB = TMOB.block(iOffB,0,this->nVB_,this->nOB_);

  RealVecMap TLinA(TExpandedA.data(),0);
  RealVecMap TLinB(TExpandedA.data(),0);

  new (&TLinA) RealVecMap(TExpandedA.data(),this->nOAVA_);
  if(doBeta) new (&TLinB) RealVecMap(TExpandedB.data(),this->nOBVB_);

  T.block(0,0,this->nOAVA_,1) = TLinA;
  if(doBeta) T.block(this->nOAVA_,0,this->nOBVB_,1) = TLinB;

} // retrieveVirOcc

template<>
void Response<double>::retrieveOccVir(RealVecMap &T, RealMatrix &TMOA,
  RealMatrix &TMOB) {

  bool doBeta = this->iPart_ != SPIN_ADAPTED && 
                this->Ref_ != SingleSlater<double>::TCS &&
                this->iClass_ != PPPA;

  auto iOffA = this->nOA_;
  auto iOffB = this->nOB_;

  RealMatrix TExpandedA, TExpandedB;
  TExpandedA = TMOA.block(0,iOffA,this->nOA_,this->nVA_);
  if(doBeta) TExpandedB = TMOB.block(0,iOffB,this->nOB_,this->nVB_);

  TExpandedA.transposeInPlace();
  if(doBeta) TExpandedB.transposeInPlace();

  RealVecMap TLinA(TExpandedA.data(),0);
  RealVecMap TLinB(TExpandedA.data(),0);

  new (&TLinA) RealVecMap(TExpandedA.data(),this->nOAVA_);
  if(doBeta) new (&TLinB) RealVecMap(TExpandedB.data(),this->nOBVB_);

  // Wierd sign that I can't consolidate
  T.block(0,0,this->nOAVA_,1) = -TLinA;
  if(doBeta) T.block(this->nOAVA_,0,this->nOBVB_,1) = -TLinB;

} // retrieveOccVir

template<>
void Response<double>::formAOTransDen(RealVecMap &T, RealMatrix &TAOA,
  RealMatrix &TAOB) {
  RealMatrix TMOA,TMOB;
  bool doBeta = this->iPart_ != SPIN_ADAPTED && 
                this->Ref_ != SingleSlater<double>::TCS &&
                this->iClass_ != PPPA;

  TMOA = RealMatrix(this->nBasis_,this->nBasis_);
  if(doBeta)
    TMOB = RealMatrix(this->nBasis_,this->nBasis_);

  if(this->iClass_ == FOPPA) {
    this->placeVirOcc(T,TMOA,TMOB);
    if(!this->doTDA_){
      RealVecMap Y(
        T.data()+this->nSingleDim_/2,this->nSingleDim_/2
      );

      this->placeOccVir(Y,TMOA,TMOB);
    }
  } else
    return;

  TAOA = (*this->singleSlater_->moA()) * TMOA * 
         this->singleSlater_->moA()->adjoint();
  if(doBeta){
    if(this->singleSlater_->isClosedShell)
      TAOB = (*this->singleSlater_->moA()) * TMOB * 
             this->singleSlater_->moA()->adjoint();
    else
      TAOB = (*this->singleSlater_->moB()) * TMOB * 
             this->singleSlater_->moB()->adjoint();
  }
}; //formAOTDen

template<>
void Response<double>::formMOTransDen(RealVecMap &T, RealMatrix &TAOA,
  RealMatrix &TAOB) {
  RealMatrix TMOA,TMOB;
  bool doBeta = this->iPart_ != SPIN_ADAPTED && 
                this->Ref_ != SingleSlater<double>::TCS &&
                this->iClass_ != PPPA;

  TMOA = RealMatrix(this->nBasis_,this->nBasis_);
  if(doBeta)
    TMOB = RealMatrix(this->nBasis_,this->nBasis_);

  TMOA = this->singleSlater_->moA()->adjoint() * TAOA * 
         (*this->singleSlater_->moA());
  if(doBeta){
    if(this->singleSlater_->isClosedShell)
      TMOB = this->singleSlater_->moA()->adjoint() * TAOB * 
             (*this->singleSlater_->moA());
    else
      TMOB = this->singleSlater_->moB()->adjoint() * TAOB * 
             (*this->singleSlater_->moB());
  }

  if(this->iClass_ == FOPPA){
    this->retrieveVirOcc(T,TMOA,TMOB);
    if(!this->doTDA_){
      RealVecMap Y(
        T.data()+this->nSingleDim_/2,this->nSingleDim_/2
      );

      this->retrieveOccVir(Y,TMOA,TMOB);
    }
  } else
    return;

}; //formMOTDen

}; // namespace ChronusQ
