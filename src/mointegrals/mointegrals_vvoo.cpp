#include <mointegrals.h>

namespace ChronusQ {
template<>
void MOIntegrals<double>::formVVOO() {
  if(this->haveMOVVOO_) return;

  bool is2C        = this->wfn_->nTCS() == 2;
  bool isOpenShell = !this->wfn_->isClosedShell;
  int NB           = this->wfn_->nBasis();
  int intDim       = this->pscf_.nOAVA() * this->pscf_.nOAVA();
  int NOA = this->wfn_->nOA();
  int NVA = this->wfn_->nVA();
  int NOB = this->wfn_->nOB();
  int NVB = this->wfn_->nVB();


  /// Evaluate VVOO (AAAA) MO Integrals

  // Allocate and zero out VVOO (AAAA) storage
  VVOOAAAA_ = this->memManager_->malloc<double>(NVA*NOA*NOA*NVA);
  std::fill_n(VVOOAAAA_,NVA*NOA*NOA*NVA,0.0);

  // Allocate and zero out Itermediates
  double * I_mnis_A   = this->memManager_->malloc<double>(NB*NB*NB*NOA);
  double * I_mnij_AA  = this->memManager_->malloc<double>(NB*NB*NOA*NOA);
  double * I_mbij_AAA = this->memManager_->malloc<double>(NB*NOA*NOA*NVA);

  std::fill_n(I_mnis_A  ,NB*NB*NB*NOA,0.0);
  std::fill_n(I_mnij_AA ,NB*NB*NOA*NOA,0.0);
  std::fill_n(I_mbij_AAA,NB*NOA*NOA*NVA,0.0);

/*
  // First Quarter (AAAA) transformation (mn|ls) -> (mn|is) (A)
  for(auto sg = 0; sg < NB; sg++)
  for(auto i = 0; i < NOA; i++)
  for(auto lm = 0; lm < NB; lm++)
  for(auto nu = 0; nu < NB; nu++)  
  for(auto mu = 0; mu < NB; mu++) {
    I_mnis_A[mu + NB*nu + NB*NB*i + NB*NB*NOA*sg] +=
      (*this->wfn_->moA())(lm,i) *
      (*this->wfn_->aointegrals()->aoERI_)(mu,nu,lm,sg);
  }

  // Second Quarter (AAAA) transformation (mn|is) (A) -> (mn|ij) (AA)
  for(auto j = 0; j < NOA; j++)
  for(auto sg = 0; sg < NB; sg++)
  for(auto i = 0; i < NOA; i++)
  for(auto nu = 0; nu < NB; nu++)  
  for(auto mu = 0; mu < NB; mu++) {
    I_mnij_AA[mu + NB*nu + NB*NB*i + NB*NB*NOA*j] +=
      (*this->wfn_->moA())(sg,j) *
      I_mnis_A[mu + NB*nu + NB*NB*i + NB*NB*NOA*sg];
  }

  // Third Quarter (AAAA) transformation (mn|ij) (AA) -> (mb|ij) (AAA)
  for(auto j = 0; j < NOA; j++)
  for(auto i = 0; i < NOA; i++)
  for(auto b = 0; b < NVA; b++)
  for(auto nu = 0; nu < NB; nu++)  
  for(auto mu = 0; mu < NB; mu++) {
    I_mbij_AAA[mu + NB*b + NVA*NB*i + NB*NVA*NOA*j] +=
      (*this->wfn_->moA())(nu,b+NOA) *
      I_mnij_AA[mu + NB*nu + NB*NB*i + NB*NB*NOA*j];
  }

  // Fourth Quarter (AAAA) transformation (mb|ij) (AAA) -> (ab|ij) (AAAA)
  for(auto j = 0; j < NOA; j++)
  for(auto i = 0; i < NOA; i++)
  for(auto b = 0; b < NVA; b++)
  for(auto a = 0; a < NVA; a++)
  for(auto mu = 0; mu < NB; mu++) {
    VVOOAAAA_[a + NVA*b + NVA*NVA*i + NVA*NVA*NOA*j] +=
      (*this->wfn_->moA())(mu,a+NOA) *
      I_mbij_AAA[mu + NB*b + NVA*NB*i + NB*NVA*NOA*j];
  }
*/
  // First Quarter (AAAA) transformation (mn|ls) -> (mn|is) (A)
  rank4w2Contract(3,this->wfn_->moA()->data(),NB,
    &this->wfn_->aointegrals()->aoERI_->storage()[0],NB,NB,NB,NB,
    I_mnis_A,NOA);

  // Second Quarter (AAAA) transformation (mn|is) (A) -> (mn|ij) (AA)
  rank4w2Contract(4,this->wfn_->moA()->data(),NB,
    I_mnis_A,NB,NB,NOA,NB,
    I_mnij_AA,NOA);

  // Third Quarter (AAAA) transformation (mn|ij) (AA) -> (mb|ij) (AAA)
  rank4w2Contract(2,this->wfn_->moA()->data() + NB*NOA,NB,
    I_mnij_AA,NB,NB,NOA,NOA,
    I_mbij_AAA,NVA);

  // Fourth Quarter (AAAA) transformation (mb|ij) (AAA) -> (ab|ij) (AAAA)
  rank4w2Contract(1,this->wfn_->moA()->data() + NB*NOA,NB,
    I_mbij_AAA,NB,NVA,NOA,NOA,
    VVOOAAAA_,NVA);

  // Free intermediate 1 and 3
  this->memManager_->free(I_mnis_A  ,NB*NB*NB*NOA);
  this->memManager_->free(I_mbij_AAA,NB*NOA*NOA*NVA);

  // If its closed shell, free intermediate 2
  if(!isOpenShell)
    this->memManager_->free(I_mnij_AA ,NB*NB*NOA*NOA);

  if(isOpenShell) {
    
    /// Evaluate VVOO (BBAA) MO Integrals

    // Allocate and Zero out VVOO (BBAA) storage
    VVOOBBAA_ = this->memManager_->malloc<double>(NVB*NOA*NOA*NVB);
    std::fill_n(VVOOBBAA_,NVB*NOA*NOA*NVB,0.0);

    // Allocate and zero out Itermediates
    double * I_mbij_BAA = this->memManager_->malloc<double>(NB*NOA*NOA*NVB);
    std::fill_n(I_mbij_BAA,NB*NOA*NOA*NVB,0.0);

/*
    // Third Quarter (BBAA) transformation (mn|ij) (AA) -> (mb|ij) (BAA)
    for(auto j = 0; j < NOA; j++)
    for(auto i = 0; i < NOA; i++)
    for(auto b = 0; b < NVB; b++)
    for(auto nu = 0; nu < NB; nu++)  
    for(auto mu = 0; mu < NB; mu++) {
      I_mbij_BAA[mu + NB*b + NVB*NB*i + NB*NVB*NOA*j] +=
        (*this->wfn_->moB())(nu,b+NOB) *
        I_mnij_AA[mu + NB*nu + NB*NB*i + NB*NB*NOA*j];
    }
 
    // Fourth Quarter (BBAA) transformation (mb|ij) (BAA) -> (ab|ij) (BBAA)
    for(auto j = 0; j < NOA; j++)
    for(auto i = 0; i < NOA; i++)
    for(auto b = 0; b < NVB; b++)
    for(auto a = 0; a < NVB; a++)
    for(auto mu = 0; mu < NB; mu++) {
      VVOOBBAA_[a + NVB*b + NVB*NVB*i + NVB*NVB*NOA*j] +=
        (*this->wfn_->moB())(mu,a+NOB) *
        I_mbij_BAA[mu + NB*b + NVB*NB*i + NB*NVB*NOA*j];
    }
*/
    // Third Quarter (BBAA) transformation (mn|ij) (AA) -> (mb|ij) (BAA)
    rank4w2Contract(2,this->wfn_->moB()->data() + NB*NOB,NB,
      I_mnij_AA,NB,NB,NOA,NOA,
      I_mbij_BAA,NVB);
    
    // Fourth Quarter (BBAA) transformation (mb|ij) (BAA) -> (ab|ij) (BBAA)
    rank4w2Contract(1,this->wfn_->moB()->data() + NB*NOB,NB,
      I_mbij_BAA,NB,NVB,NOA,NOA,
      VVOOBBAA_,NVB);


    // Free intermediates 2 and 3
    this->memManager_->free(I_mnij_AA ,NB*NB*NOA*NOA);
    this->memManager_->free(I_mbij_BAA,NB*NOA*NOA*NVB);

    /// Evaluate VVOO (BBBB) MO Integrals

    // Allocate and zero out VVOO (BBBB) storage
    VVOOBBBB_ = this->memManager_->malloc<double>(NVB*NOB*NOB*NVB);
    std::fill_n(VVOOBBBB_,NVB*NOB*NOB*NVB,0.0);

    // Allocate and zero out Itermediates
    double * I_mnis_B   = this->memManager_->malloc<double>(NB*NB*NB*NOB);
    double * I_mnij_BB  = this->memManager_->malloc<double>(NB*NB*NOB*NOB);
    double * I_mbij_BBB = this->memManager_->malloc<double>(NB*NOB*NOB*NVB);

    std::fill_n(I_mnis_B  ,NB*NB*NB*NOB,0.0);
    std::fill_n(I_mnij_BB ,NB*NB*NOB*NOB,0.0);
    std::fill_n(I_mbij_BBB,NB*NOB*NOB*NVB,0.0);

/*
    for(auto sg = 0; sg < NB; sg++)
    for(auto i = 0; i < NOB; i++)
    for(auto lm = 0; lm < NB; lm++)
    for(auto nu = 0; nu < NB; nu++)  
    for(auto mu = 0; mu < NB; mu++) {
      I_mnis_B[mu + NB*nu + NB*NB*i + NB*NB*NOB*sg] +=
        (*this->wfn_->moB())(lm,i) *
        (*this->wfn_->aointegrals()->aoERI_)(mu,nu,lm,sg);
    }

    for(auto j = 0; j < NOB; j++)
    for(auto sg = 0; sg < NB; sg++)
    for(auto i = 0; i < NOB; i++)
    for(auto nu = 0; nu < NB; nu++)  
    for(auto mu = 0; mu < NB; mu++) {
      I_mnij_BB[mu + NB*nu + NB*NB*i + NB*NB*NOB*j] +=
        (*this->wfn_->moB())(sg,j) *
        I_mnis_B[mu + NB*nu + NB*NB*i + NB*NB*NOB*sg];
    }

    for(auto j = 0; j < NOB; j++)
    for(auto i = 0; i < NOB; i++)
    for(auto b = 0; b < NVB; b++)
    for(auto nu = 0; nu < NB; nu++)  
    for(auto mu = 0; mu < NB; mu++) {
      I_mbij_BBB[mu + NB*b + NVB*NB*i + NB*NVB*NOB*j] +=
        (*this->wfn_->moB())(nu,b+NOB) *
        I_mnij_BB[mu + NB*nu + NB*NB*i + NB*NB*NOB*j];
    }

    for(auto j = 0; j < NOB; j++)
    for(auto i = 0; i < NOB; i++)
    for(auto b = 0; b < NVB; b++)
    for(auto a = 0; a < NVB; a++)
    for(auto mu = 0; mu < NB; mu++) {
      VVOOBBBB_[a + NVB*b + NVB*NVB*i + NVB*NVB*NOB*j] +=
        (*this->wfn_->moB())(mu,a+NOB) *
        I_mbij_BBB[mu + NB*b + NVB*NB*i + NB*NVB*NOB*j];
    }
*/
    // First Quarter (BBBB) transformation (mn|ls) -> (mn|is) (B)
    rank4w2Contract(3,this->wfn_->moB()->data(),NB,
      &this->wfn_->aointegrals()->aoERI_->storage()[0],NB,NB,NB,NB,
      I_mnis_B,NOB);

    // Second Quarter (BBBB) transformation (mn|is) (B) -> (mn|ij) (BB)
    rank4w2Contract(4,this->wfn_->moB()->data(),NB,
      I_mnis_B,NB,NB,NOB,NB,
      I_mnij_BB,NOB);

    // Third Quarter (BBBB) transformation (mn|ij) (BB) -> (mb|ij) (BBB)
    rank4w2Contract(2,this->wfn_->moB()->data() + NB*NOB,NB,
      I_mnij_BB,NB,NB,NOB,NOB,
      I_mbij_BBB,NVB);

    // Fourth Quarter (BBBB) transformation (mb|ij) (BBB) -> (ab|ij) (BBBB)
    rank4w2Contract(1,this->wfn_->moB()->data() + NB*NOB,NB,
      I_mbij_BBB,NB,NVB,NOB,NOB,
      VVOOBBBB_,NVB);

    this->memManager_->free(I_mnis_B  ,NB*NB*NB*NOB);
    this->memManager_->free(I_mbij_BBB,NB*NOB*NOB*NVB);

    VVOOAABB_ = this->memManager_->malloc<double>(NVA*NOB*NOB*NVA);
    std::fill_n(VVOOAABB_,NVA*NOB*NOB*NVA,0.0);

    double * I_mbij_ABB = this->memManager_->malloc<double>(NB*NOB*NOB*NVA);
    std::fill_n(I_mbij_ABB,NOB*NOB*NVA*NB,0.0);

/*
    for(auto j = 0; j < NOB; j++)
    for(auto i = 0; i < NOB; i++)
    for(auto b = 0; b < NVA; b++)
    for(auto nu = 0; nu < NB; nu++)  
    for(auto mu = 0; mu < NB; mu++) {
      I_mbij_ABB[mu + NB*b + NVA*NB*i + NB*NVA*NOB*j] +=
        (*this->wfn_->moA())(nu,b+NOA) *
        I_mnij_BB[mu + NB*nu + NB*NB*i + NB*NB*NOB*j];
    }

    for(auto j = 0; j < NOB; j++)
    for(auto i = 0; i < NOB; i++)
    for(auto b = 0; b < NVA; b++)
    for(auto a = 0; a < NVA; a++)
    for(auto mu = 0; mu < NB; mu++) {
      VVOOAABB_[a + NVA*b + NVA*NVA*i + NVA*NVA*NOB*j] +=
        (*this->wfn_->moA())(mu,a+NOA) *
        I_mbij_ABB[mu + NB*b + NVA*NB*i + NB*NVA*NOB*j];
    }
*/

    // Third Quarter (AABB) transformation (mn|ij) (BB) -> (mb|ij) (ABB)
    rank4w2Contract(2,this->wfn_->moA()->data() + NB*NOA,NB,
      I_mnij_BB,NB,NB,NOB,NOB,
      I_mbij_ABB,NVA);
    
    // Fourth Quarter (AABB) transformation (mb|ij) (ABB) -> (ab|ij) (AABB)
    rank4w2Contract(1,this->wfn_->moA()->data() + NB*NOA,NB,
      I_mbij_ABB,NB,NVA,NOB,NOB,
      VVOOAABB_,NVB);

    this->memManager_->free(I_mnij_BB ,NB*NB*NOB*NOB);
    this->memManager_->free(I_mbij_ABB,NB*NOB*NOB*NVA);
  }

  // Realize that spin blocks are equivalent for RHF
  if(!isOpenShell) {
    VVOOAABB_ = VVOOAAAA_;
    VVOOBBAA_ = VVOOAAAA_;
    VVOOBBBB_ = VVOOAAAA_;
  }

  this->haveMOVVOO_ = true;
} // formVVOO

template<>
void MOIntegrals<double>::formFullVVOO(){
  this->formVVOO();
  int NOA = this->wfn_->nOA();
  int NVA = this->wfn_->nVA();
  int NOB = this->wfn_->nOB();
  int NVB = this->wfn_->nVB();
  int NO = this->wfn_->nO();
  int NV = this->wfn_->nV();

  this->VVOO_ = this->memManager_->malloc<double>(NV*NO*NV*NO);
  std::fill_n(this->VVOO_,NO*NV*NO*NV,0.0);

  for(auto J = 0; J < NO; J+=2)
  for(auto I = 0; I < NO; I+=2)
  for(auto B = 0; B < NV; B+=2)
  for(auto A = 0; A < NV; A+=2){
    int a = A / 2;
    int i = I / 2;
    int b = B / 2;
    int j = J / 2;
    this->VVOO_[A + B*NV + I*NV*NV + J*NO*NV*NV] = 
      this->VVOOAAAA_[a + b*NVA + i*NVA*NVA + j*NOA*NVA*NVA];

    this->VVOO_[(A+1) + (B+1)*NV + I*NV*NV + J*NO*NV*NV] = 
      this->VVOOBBAA_[a + b*NVB + i*NVB*NVB + j*NOA*NVB*NVB];

    this->VVOO_[A + B*NV + (I+1)*NV*NV + (J+1)*NO*NV*NV] = 
      this->VVOOAABB_[a + b*NVA + i*NVA*NVA + j*NOB*NVA*NVA];

    this->VVOO_[(A+1) + (B+1)*NV + (I+1)*NV*NV + (J+1)*NO*NV*NV] = 
      this->VVOOBBBB_[a + b*NVB + i*NVB*NVB + j*NOB*NVB*NVB];
  }
/*
  for(auto j = 0; j < NOA; j++)
  for(auto i = 0; i < NOA; i++)
  for(auto b = 0; b < NVA; b++)
  for(auto a = 0; a < NVA; a++){
    this->VVOO_[a + b*NV + i*NV*NV + j*NO*NV*NV] = 
      this->VVOOAAAA_[a + b*NVA + i*NVA*NVA + j*NOA*NVA*NVA];
  }
  for(auto j = 0; j < NOB; j++)
  for(auto i = 0; i < NOB; i++)
  for(auto b = 0; b < NVA; b++)
  for(auto a = 0; a < NVA; a++){
    this->VVOO_[a + b*NV + (i+NOA)*NV*NV + (j+NOA)*NO*NV*NV] = 
      this->VVOOAABB_[a + b*NVA + i*NVA*NVA + j*NOB*NVA*NVA];
  }
  for(auto j = 0; j < NOA; j++)
  for(auto i = 0; i < NOA; i++)
  for(auto b = 0; b < NVB; b++)
  for(auto a = 0; a < NVB; a++){
    this->VVOO_[(a+NVA) + (b+NVA)*NV + i*NV*NV + j*NO*NV*NV] = 
      this->VVOOBBAA_[a + b*NVB + i*NVB*NVB + i*NOA*NVB*NVB];
  }
  for(auto j = 0; j < NOB; j++)
  for(auto i = 0; i < NOB; i++)
  for(auto b = 0; b < NVB; b++)
  for(auto a = 0; a < NVB; a++){
    this->VVOO_[(a+NVA) + (b+NVA)*NV + (i+NOA)*NV*NV + (j+NOA)*NO*NV*NV] = 
      this->VVOOBBBB_[a + i*NVB + b*NVB*NVB + j*NOB*NVB*NVB];
  }
*/
} // formFullVVOO

};
