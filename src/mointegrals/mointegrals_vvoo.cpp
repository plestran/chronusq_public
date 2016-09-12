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


  /** Evaluate VVOO (AAAA) MO Integrals **/

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

  // Free intermediate 1 and 3
  this->memManager_->free(I_mnis_A  ,NB*NB*NB*NOA);
  this->memManager_->free(I_mbij_AAA,NB*NOA*NOA*NVA);

  // If its closed shell, free intermediate 2
  if(!isOpenShell)
    this->memManager_->free(I_mnij_AA ,NB*NB*NOA*NOA);

  if(isOpenShell) {
    
    /** Evaluate VVOO (BBAA) MO Integrals **/

    // Allocate and Zero out VVOO (BBAA) storage
    VVOOBBAA_ = this->memManager_->malloc<double>(NVB*NOA*NOA*NVB);
    std::fill_n(VVOOAAAA_,NVB*NOA*NOA*NVB,0.0);

    // Allocate and zero out Itermediates
    double * I_mbij_BAA = this->memManager_->malloc<double>(NB*NOA*NOA*NVB);
    std::fill_n(I_mbij_BAA,NB*NOA*NOA*NVB,0.0);

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


    // Free intermediates 2 and 3
    this->memManager_->free(I_mnij_AA ,NB*NB*NOA*NOA);
    this->memManager_->free(I_mbij_BAA,NB*NOA*NOA*NVB);

    /** Evaluate VVOO (BBBB) MO Integrals **/

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

    this->memManager_->free(I_mnis_B  ,NB*NB*NB*NOB);
    this->memManager_->free(I_mbij_BBB,NB*NOB*NOB*NVB);

    VVOOAABB_ = this->memManager_->malloc<double>(NVA*NOB*NOB*NVA);
    std::fill_n(VVOOAABB_,NVA*NOB*NOB*NVA,0.0);

    double * I_mbij_ABB = this->memManager_->malloc<double>(NB*NOB*NOB*NVA);
    std::fill_n(I_mbij_ABB,NOB*NOB*NVA*NB,0.0);

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


    this->memManager_->free(I_mnij_BB ,NB*NB*NOB*NOB);
    this->memManager_->free(I_mbij_ABB,NB*NOB*NOB*NVA);
  }


  this->haveMOVVOO_ = true;
} // formVVOO

template<>
void MOIntegrals<double>::formFullVVOO(){
  this->formVVOO();
  int NOA = this->wfn_->nOA();
  int NVA = this->wfn_->nVA();
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
      this->VVOOAAAA_[a + b*NVA + i*NVA*NVA + j*NOA*NVA*NVA];

    this->VVOO_[A + B*NV + (I+1)*NV*NV + (J+1)*NO*NV*NV] = 
      this->VVOOAAAA_[a + b*NVA + i*NVA*NVA + j*NOA*NVA*NVA];

    this->VVOO_[(A+1) + (B+1)*NV + (I+1)*NV*NV + (J+1)*NO*NV*NV] = 
      this->VVOOAAAA_[a + b*NVA + i*NVA*NVA + j*NOA*NVA*NVA];
  }
} // formFullVVOO

};
