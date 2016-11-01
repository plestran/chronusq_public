#include <mointegrals.h>

namespace ChronusQ {
template<>
void MOIntegrals<double>::formVOOO(){
  if(this->haveMOVOOO_) return;

  bool is2C        = this->wfn_->nTCS() == 2;
  bool isOpenShell = !this->wfn_->isClosedShell;
  int NB           = this->wfn_->nBasis();
  int intDim       = this->pscf_.nOAVA() * this->pscf_.nOAVA();
  int NOA = this->wfn_->nOA();
  int NVA = this->wfn_->nVA();
  int NOB = this->wfn_->nOB();
  int NVB = this->wfn_->nVB();


  /// Evaluate VOOO (AAAA) MO Integrals

  // Allocate and zero out VOOO (AAAA) storage
  VOOOAAAA_ = this->memManager_->malloc<double>(NVA*NOA*NOA*NOA);
  std::fill_n(VOOOAAAA_,NVA*NOA*NOA*NOA,0.0);

  // Allocate and zero out Itermediates
  double * I_1_A   = this->memManager_->malloc<double>(NB*NB*NB*NOA);
  double * I_2_AA  = this->memManager_->malloc<double>(NB*NB*NOA*NOA);
  double * I_3_AAA = this->memManager_->malloc<double>(NB*NOA*NOA*NOA);

  std::fill_n(I_1_A  ,NB*NB*NB*NOA,0.0);
  std::fill_n(I_2_AA ,NB*NB*NOA*NOA,0.0);
  std::fill_n(I_3_AAA,NB*NOA*NOA*NOA,0.0);

  // First Quarter (AAAA) transformation (mn|ls) -> (mn|lk) (A)
  rank4w2Contract(4,this->wfn_->moA()->data(),NB,
    &this->wfn_->aointegrals()->aoERI_->storage()[0],NB,NB,NB,NB,
    I_1_A,NOA);

  // Second Quarter (AAAA) transformation (mn|lk) (A) -> (mn|jk) (AA)
  rank4w2Contract(3,this->wfn_->moA()->data(),NB,  
    I_1_A,NB,NB,NB,NOA,
    I_2_AA,NOA);

  // Third Quarter (AAAA) transformation (mn|jk) (AA) -> (mi|jk) (AAA)
  rank4w2Contract(2,this->wfn_->moA()->data(),NB,  
    I_2_AA,NB,NB,NOA,NOA,
    I_3_AAA,NOA);

  // Fourth Quarter (AAAA) transformation (mi|jk) (AAA) -> (ai|jk) (AAAA)
  rank4w2Contract(1,this->wfn_->moA()->data() + NOA*NB,NB,  
    I_3_AAA,NB,NOA,NOA,NOA,
    VOOOAAAA_,NVA);

  this->memManager_->free(I_3_AAA,NB*NOA*NOA*NOA);
  this->memManager_->free(I_1_A  ,NB*NB*NB*NOA);

  if(isOpenShell) {
    // Allocate and zero out VOOO (BBAA) storage
    VOOOBBAA_ = this->memManager_->malloc<double>(NVB*NOB*NOA*NOA);
    std::fill_n(VOOOBBAA_,NVB*NOB*NOA*NOA,0.0);

    // Allocate and zero out Itermediates
    double * I_3_BAA = this->memManager_->malloc<double>(NB*NOB*NOA*NOA);

    std::fill_n(I_3_BAA,NB*NOB*NOA*NOA,0.0);

    // Third Quarter (AAAA) transformation (mn|jk) (AA) -> (mi|jk) (BAA)
    rank4w2Contract(2,this->wfn_->moB()->data(),NB,  
      I_2_AA,NB,NB,NOA,NOA,
      I_3_BAA,NOB);
 
    // Fourth Quarter (AAAA) transformation (mi|jk) (BAA) -> (ai|jk) (BBAA)
    rank4w2Contract(1,this->wfn_->moB()->data() + NOB*NB,NB,  
      I_3_BAA,NB,NOB,NOA,NOA,
      VOOOBBAA_,NVB);

    this->memManager_->free(I_3_BAA,NB*NOB*NOA*NOA);

  }

  this->memManager_->free(I_2_AA ,NB*NB*NOA*NOA);

  if(isOpenShell) {
    /// Evaluate VOOO (BBBB) MO Integrals
 
    // Allocate and zero out VOOO (BBBB) storage
    VOOOBBBB_ = this->memManager_->malloc<double>(NVB*NOB*NOB*NOB);
    std::fill_n(VOOOBBBB_,NVB*NOB*NOB*NOB,0.0);
 
    // Allocate and zero out Itermediates
    double * I_1_B   = this->memManager_->malloc<double>(NB*NB*NB*NOB);
    double * I_2_BB  = this->memManager_->malloc<double>(NB*NB*NOB*NOB);
    double * I_3_BBB = this->memManager_->malloc<double>(NB*NOB*NOB*NOB);
 
    std::fill_n(I_1_B  ,NB*NB*NB*NOB,0.0);
    std::fill_n(I_2_BB ,NB*NB*NOB*NOB,0.0);
    std::fill_n(I_3_BBB,NB*NOB*NOB*NOB,0.0);
 
    // First Quarter (BBBB) transformation (mn|ls) -> (mn|lk) (B)
    rank4w2Contract(4,this->wfn_->moB()->data(),NB,
      &this->wfn_->aointegrals()->aoERI_->storage()[0],NB,NB,NB,NB,
      I_1_B,NOB);
 
    // Second Quarter (BBBB) transformation (mn|lk) (B) -> (mn|jk) (BB)
    rank4w2Contract(3,this->wfn_->moB()->data(),NB,  
      I_1_B,NB,NB,NB,NOB,
      I_2_BB,NOB);
 
    // Third Quarter (BBBB) transformation (mn|jk) (BB) -> (mi|jk) (BBB)
    rank4w2Contract(2,this->wfn_->moB()->data(),NB,  
      I_2_BB,NB,NB,NOB,NOB,
      I_3_BBB,NOB);
 
    // Fourth Quarter (BBBB) transformation (mi|jk) (BBB) -> (ai|jk) (BBBB)
    rank4w2Contract(1,this->wfn_->moB()->data() + NOB*NB,NB,  
      I_3_BBB,NB,NOB,NOB,NOB,
      VOOOBBBB_,NVB);
 
    this->memManager_->free(I_3_BBB,NB*NOB*NOB*NOB);
    this->memManager_->free(I_1_B  ,NB*NB*NB*NOB);


    // Allocate and zero out VOOO (AABB) storage
    VOOOAABB_ = this->memManager_->malloc<double>(NVA*NOA*NOB*NOB);
    std::fill_n(VOOOAABB_,NVA*NOA*NOB*NOB,0.0);

    // Allocate and zero out Itermediates
    double * I_3_ABB = this->memManager_->malloc<double>(NB*NOA*NOB*NOB);

    std::fill_n(I_3_ABB,NB*NOA*NOB*NOB,0.0);

    // Third Quarter (AAAA) transformation (mn|jk) (BB) -> (mi|jk) (ABB)
    rank4w2Contract(2,this->wfn_->moA()->data(),NB,  
      I_2_AA,NB,NB,NOB,NOB,
      I_3_ABB,NOA);
 
    // Fourth Quarter (AAAA) transformation (mi|jk) (ABB) -> (ai|jk) (AABB)
    rank4w2Contract(1,this->wfn_->moA()->data() + NOA*NB,NB,  
      I_3_ABB,NB,NOA,NOB,NOB,
      VOOOAABB_,NVA);

    this->memManager_->free(I_3_ABB,NB*NOA*NOB*NOB);
    this->memManager_->free(I_2_BB ,NB*NB*NOB*NOB);
  }
}

template<>
void MOIntegrals<double>::formFullVOOO(){
  this->formVOOO();
  int NOA = this->wfn_->nOA();
  int NVA = this->wfn_->nVA();
  int NOB = this->wfn_->nOB();
  int NVB = this->wfn_->nVB();
  int NO = this->wfn_->nO();
  int NV = this->wfn_->nV();

  this->VOOO_ = this->memManager_->malloc<double>(NV*NO*NO*NO);
  std::fill_n(this->VOOO_,NO*NO*NO*NV,0.0);

  for(auto K = 0; K < NO; K+=2)
  for(auto J = 0; J < NO; J+=2)
  for(auto I = 0; I < NO; I+=2)
  for(auto A = 0; A < NV; A+=2){
    int a = A / 2;
    int i = I / 2;
    int j = J / 2;
    int k = K / 2;
    this->VOOO_[A + I*NV + J*NV*NV + K*NO*NV*NV] = 
      this->VOOOAAAA_[a + i*NVA + j*NVA*NOA + k*NVA*NOA*NOA];

    this->VOOO_[(A+1) + (I+1)*NV + J*NV*NV + K*NO*NV*NV] = 
      this->VOOOBBAA_[a + i*NVB + j*NVB*NOB + k*NVB*NOB*NOA];

    this->VOOO_[A + I*NV + (J+1)*NV*NV + (K+1)*NO*NV*NV] = 
      this->VOOOAABB_[a + j*NVA + j*NVA*NOA + k*NVA*NOA*NOB];

    this->VOOO_[(A+1) + (I+1)*NV + (J+1)*NV*NV + (K+1)*NO*NV*NV] = 
      this->VOOOBBBB_[a + k*NVB + j*NVB*NOB + k*NVB*NOB*NOB];
  }
} // formFullVOOO
};
