#include <mointegrals.h>

namespace ChronusQ {
template<>
void MOIntegrals<double>::formVVVV() {
  if(this->haveMOVVVV_) return;

  bool is2C        = this->wfn_->nTCS() == 2;
  bool isOpenShell = !this->wfn_->isClosedShell;
  int NB           = this->wfn_->nBasis();
  int intDim       = this->pscf_.nVAVA() * this->pscf_.nVAVA();
  int NVA = this->wfn_->nVA();
  int NOA = this->wfn_->nOA();
  int NVB = this->wfn_->nVB();
  int NOB = this->wfn_->nOB();


  /// Evaluate VVVV (AAAA) MO Integrals

  // Allocate and zero out VVVV (AAAA) storage
  VVVVAAAA_ = this->memManager_->malloc<double>(NVA*NVA*NVA*NVA);
  std::fill_n(VVVVAAAA_,NVA*NVA*NVA*NVA,0.0);

  // Allocate and zero out Itermediates
  double * I_anls_A   = this->memManager_->malloc<double>(NB*NB*NB*NVA);
  double * I_abls_AA  = this->memManager_->malloc<double>(NB*NB*NVA*NVA);
  double * I_abcs_AAA = this->memManager_->malloc<double>(NB*NVA*NVA*NVA);

  std::fill_n(I_anls_A  ,NB*NB*NB*NVA,0.0);
  std::fill_n(I_abls_AA ,NB*NB*NVA*NVA,0.0);
  std::fill_n(I_abcs_AAA,NB*NVA*NVA*NVA,0.0);

  // First Quarter (AAAA) transformation (mn|ls) -> (an|ls) (A)
  rank4w2Contract(1,this->wfn_->moA()->data() + NOA*NB,NB,
    &this->wfn_->aointegrals()->aoERI_->storage()[0],NB,NB,NB,NB,
    I_anls_A,NVA);

  // Second Quarter (AAAA) transformation (an|ls) (A) -> (ab|ls) (AA)
  rank4w2Contract(2,this->wfn_->moA()->data() + NOA*NB,NB,  
    I_anls_A,NVA,NB,NB,NB,
    I_abls_AA,NVA);

  // Third Quarter (AAAA) transformation (ab|ls) (AA) -> (ab|cs) (AAA)
  rank4w2Contract(3,this->wfn_->moA()->data() + NOA*NB,NB,  
    I_abls_AA,NVA,NVA,NB,NB,
    I_abcs_AAA,NVA);

  // Fourth Quarter (AAAA) transformation (ab|cs) (AAA) -> (ab|cd) (AAAA)
  rank4w2Contract(4,this->wfn_->moA()->data() + NOA*NB,NB,  
    I_abcs_AAA,NVA,NVA,NVA,NB,
    VVVVAAAA_,NVA);

  this->memManager_->free(I_abcs_AAA,NB*NVA*NVA*NVA);
  this->memManager_->free(I_anls_A  ,NB*NB*NB*NVA);

  if(isOpenShell) {
    /// Evaluate VVVV (AABB) MO Integrals

    // Allocate and zero out VVVV (AABB) storage
    VVVVAABB_ = this->memManager_->malloc<double>(NVA*NVA*NVB*NVB);
    std::fill_n(VVVVAABB_,NVA*NVA*NVB*NVB,0.0);

    // Allocate and zero out Itermediates
    double * I_abcs_AAB = this->memManager_->malloc<double>(NB*NVA*NVA*NVB);

    std::fill_n(I_abcs_AAB,NB*NVA*NVA*NVB,0.0);

    // Third Quarter (AABB) transformation (ab|ls) (AA) -> (ab|cs) (AAB)
    rank4w2Contract(3,this->wfn_->moB()->data() + NOB*NB,NB,  
      I_abls_AA,NVA,NVA,NB,NB,
      I_abcs_AAB,NVB);
 
    // Fourth Quarter (AABB) transformation (ab|cs) (AAB) -> (ab|cd) (AABB)
    rank4w2Contract(4,this->wfn_->moB()->data() + NOB*NB,NB,  
      I_abcs_AAB,NVA,NVA,NVB,NB,
      VVVVAABB_,NVB);
    
    this->memManager_->free(I_abcs_AAB,NB*NVA*NVA*NVB);
  }

  this->memManager_->free(I_abls_AA ,NB*NB*NVA*NVA);

  if(isOpenShell) {
    /// Evaluate VVVV (BBBB) MO Integrals
 
    // Allocate and zero out VVVV (BBBB) storage
    VVVVBBBB_ = this->memManager_->malloc<double>(NVB*NVB*NVB*NVB);
    std::fill_n(VVVVBBBB_,NVB*NVB*NVB*NVB,0.0);
 
    // Allocate and zero out Itermediates
    double * I_anls_B   = this->memManager_->malloc<double>(NB*NB*NB*NVB);
    double * I_abls_BB  = this->memManager_->malloc<double>(NB*NB*NVB*NVB);
    double * I_abcs_BBB = this->memManager_->malloc<double>(NB*NVB*NVB*NVB);
 
    std::fill_n(I_anls_B  ,NB*NB*NB*NVB,0.0);
    std::fill_n(I_abls_BB ,NB*NB*NVB*NVB,0.0);
    std::fill_n(I_abcs_BBB,NB*NVB*NVB*NVB,0.0);
 
    // First Quarter (BBBB) transformation (mn|ls) -> (an|ls) (B)
    rank4w2Contract(1,this->wfn_->moB()->data() + NOB*NB,NB,
      &this->wfn_->aointegrals()->aoERI_->storage()[0],NB,NB,NB,NB,
      I_anls_B,NVB);
 
    // Second Quarter (BBBB) transformation (an|ls) (B) -> (ab|ls) (BB)
    rank4w2Contract(2,this->wfn_->moB()->data() + NOB*NB,NB,  
      I_anls_B,NVB,NB,NB,NB,
      I_abls_BB,NVB);
 
    // Third Quarter (BBBB) transformation (ab|ls) (BB) -> (ab|cs) (BBB)
    rank4w2Contract(3,this->wfn_->moB()->data() + NOB*NB,NB,  
      I_abls_BB,NVB,NVB,NB,NB,
      I_abcs_BBB,NVB);
 
    // Fourth Quarter (BBBB) transformation (ab|cs) (BBB) -> (ab|cd) (BBBB)
    rank4w2Contract(4,this->wfn_->moB()->data() + NOB*NB,NB,  
      I_abcs_BBB,NVB,NVB,NVB,NB,
      VVVVBBBB_,NVB);
 
    this->memManager_->free(I_abcs_BBB,NB*NVB*NVB*NVB);
    this->memManager_->free(I_abls_BB ,NB*NB*NVB*NVB);
    this->memManager_->free(I_anls_B  ,NB*NB*NB*NVB);
  }

  this->haveMOVVVV_ = true;
} // formVVVV

template<>
void MOIntegrals<double>::formFullVVVV(){
  this->formVVVV();
  int NOA = this->wfn_->nOA();
  int NVA = this->wfn_->nVA();
  int NOB = this->wfn_->nOB();
  int NVB = this->wfn_->nVB();
  int NO = this->wfn_->nO();
  int NV = this->wfn_->nV();

  this->VVVV_ = this->memManager_->malloc<double>(NV*NV*NV*NV);
  std::fill_n(this->VVVV_,NV*NV*NV*NV,0.0);

  for(auto D = 0; D < NV; D+=2)
  for(auto C = 0; C < NV; C+=2)
  for(auto B = 0; B < NV; B+=2)
  for(auto A = 0; A < NV; A+=2){
    int a = A / 2;
    int b = B / 2;
    int c = C / 2;
    int d = D / 2;

    this->VVVV_[A + B*NV + C*NV*NV + D*NV*NV*NV] = 
      this->VVVVAAAA_[a + b*NVA + c*NVA*NVA + d*NVA*NVA*NVA];

    this->VVVV_[(A+1) + (B+1)*NV + C*NV*NV + D*NV*NV*NV] = 
      this->VVVVAABB_[c + d*NVA + a*NVA*NVA + b*NVA*NVA*NVB];

    this->VVVV_[A + B*NV + (C+1)*NV*NV + (D+1)*NV*NV*NV] = 
      this->VVVVAABB_[a + b*NVA + c*NVA*NVA + d*NVA*NVA*NVB];

    this->VVVV_[(A+1) + (B+1)*NV + (C+1)*NV*NV + (D+1)*NV*NV*NV] = 
      this->VVVVBBBB_[a + b*NVB + c*NVB*NVB + d*NVB*NVB*NVB];
  }
} // formFullVVVV

};
