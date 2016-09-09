#include <mointegrals.h>

namespace ChronusQ {
template<>
void MOIntegrals<double>::formVOVO() {
  if(this->haveMOVOVO_) return;

  bool is2C        = this->wfn_->nTCS() == 2;
  bool isOpenShell = !this->wfn_->isClosedShell;
  int NB           = this->wfn_->nBasis();
  int intDim       = this->pscf_.nOAVA() * this->pscf_.nOAVA();
  int NOA = this->wfn_->nOA();
  int NVA = this->wfn_->nVA();
  int NOB = this->wfn_->nOB();
  int NVB = this->wfn_->nVB();


  VOVOAAAA_ = this->memManager_->malloc<double>(NOA*NVA*NOA*NVA);
  std::fill_n(VOVOAAAA_,NVA*NOA*NOA*NVA,0.0);

  double * I_mils_A   = this->memManager_->malloc<double>(NB*NB*NB*NOA);
  double * I_milj_AA  = this->memManager_->malloc<double>(NB*NB*NOA*NOA);
  double * I_mibj_AAA = this->memManager_->malloc<double>(NB*NOA*NOA*NVA);

  std::fill_n(I_mils_A  ,NB*NB*NB*NOA,0.0);
  std::fill_n(I_milj_AA ,NB*NB*NOA*NOA,0.0);
  std::fill_n(I_mibj_AAA,NB*NOA*NOA*NVA,0.0);

  for(auto sg = 0; sg < NB; sg++)
  for(auto lm = 0; lm < NB; lm++)
  for(auto i = 0; i < NOA; i++)
  for(auto nu = 0; nu < NB; nu++)  
  for(auto mu = 0; mu < NB; mu++) {
    I_mils_A[mu + NB*i + NB*NOA*lm + NB*NB*NOA*sg] +=
      (*this->wfn_->moA())(nu,i) *
      (*this->wfn_->aointegrals()->aoERI_)(mu,nu,lm,sg);
  }

  for(auto j = 0; j < NOA; j++)
  for(auto sg = 0; sg < NB; sg++)
  for(auto lm = 0; lm < NB; lm++)
  for(auto i = 0; i < NOA; i++)
  for(auto mu = 0; mu < NB; mu++) {
    I_milj_AA[mu + NB*i + NB*NOA*lm + NB*NB*NOA*j] +=
      (*this->wfn_->moA())(sg,j) *
      I_mils_A[mu + NB*i + NB*NOA*lm + NB*NB*NOA*sg];
  }

  for(auto j = 0; j < NOA; j++)
  for(auto b = 0; b < NVA; b++)
  for(auto lm = 0; lm < NB; lm++)
  for(auto i = 0; i < NOA; i++)
  for(auto mu = 0; mu < NB; mu++) {
    I_mibj_AAA[mu + NB*i + NB*NOA*b + NB*NVA*NOA*j] +=
      (*this->wfn_->moA())(lm,b+NOA) *
      I_milj_AA[mu + NB*i + NB*NOA*lm + NB*NB*NOA*j];
  }

  for(auto j = 0; j < NOA; j++)
  for(auto b = 0; b < NVA; b++)
  for(auto i = 0; i < NOA; i++)
  for(auto a = 0; a < NVA; a++)
  for(auto mu = 0; mu < NB; mu++) {
    VOVOAAAA_[a + NVA*i + NVA*NOA*b + NVA*NVA*NOA*j] +=
      (*this->wfn_->moA())(mu,a+NOA) *
      I_mibj_AAA[mu + NB*i + NB*NOA*b + NB*NVA*NOA*j];
  }


  // Deallocate the unused intermediates
  this->memManager_->free(I_mibj_AAA,NB*NOA*NOA*NVA);
  this->memManager_->free(I_milj_AA ,NB*NB*NOA*NOA);

  if(isOpenShell) {
    VOVOAABB_ = this->memManager_->malloc<double>(NOA*NVA*NOB*NVB);
    std::fill_n(VOVOAABB_,NOA*NVA*NOB*NVB,0.0);

    double * I_milj_AB  = this->memManager_->malloc<double>(NB*NB*NOA*NOB);
    double * I_mibj_ABB = this->memManager_->malloc<double>(NB*NOA*NOB*NVB);

    std::fill_n(I_milj_AB ,NB*NB*NOA*NOB,0.0);
    std::fill_n(I_mibj_ABB,NB*NOA*NOB*NVB,0.0);

    for(auto j = 0; j < NOB; j++)
    for(auto sg = 0; sg < NB; sg++)
    for(auto lm = 0; lm < NB; lm++)
    for(auto i = 0; i < NOA; i++)
    for(auto mu = 0; mu < NB; mu++) {
      I_milj_AB[mu + NB*i + NB*NOA*lm + NB*NB*NOA*j] +=
        (*this->wfn_->moB())(sg,j) *
        I_mils_A[mu + NB*i + NB*NOA*lm + NB*NB*NOA*sg];
    }
 
    for(auto j = 0; j < NOB; j++)
    for(auto b = 0; b < NVB; b++)
    for(auto lm = 0; lm < NB; lm++)
    for(auto i = 0; i < NOA; i++)
    for(auto mu = 0; mu < NB; mu++) {
      I_mibj_ABB[mu + NB*i + NB*NOA*b + NB*NVB*NOA*j] +=
        (*this->wfn_->moB())(lm,b+NOB) *
        I_milj_AB[mu + NB*i + NB*NOA*lm + NB*NB*NOA*j];
    }

    for(auto j = 0; j < NOB; j++)
    for(auto b = 0; b < NVB; b++)
    for(auto i = 0; i < NOA; i++)
    for(auto a = 0; a < NVA; a++)
    for(auto mu = 0; mu < NB; mu++) {
      VOVOAABB_[a + NVA*i + NVA*NOA*b + NVA*NVB*NOA*j] +=
        (*this->wfn_->moA())(mu,a+NOA) *
        I_mibj_ABB[mu + NB*i + NB*NOA*b + NB*NVB*NOA*j];
    }

    // Deallocate the unused intermediates
    this->memManager_->free(I_milj_AB ,NB*NB*NOA*NOB);
    this->memManager_->free(I_mibj_ABB,NB*NOA*NOB*NVB);
  }

  // Deallocate the unused intermediates
  this->memManager_->free(I_mils_A  ,NB*NB*NB*NOA);

  if(isOpenShell) {
    VOVOBBBB_ = this->memManager_->malloc<double>(NOB*NVB*NOB*NVB);
    std::fill_n(VOVOBBBB_,NVB*NOB*NOB*NVB,0.0);

    double * I_mils_B   = this->memManager_->malloc<double>(NB*NB*NB*NOB);
    double * I_milj_BB  = this->memManager_->malloc<double>(NB*NB*NOB*NOB);
    double * I_mibj_BBB = this->memManager_->malloc<double>(NB*NOB*NOB*NVB);

    std::fill_n(I_mils_B  ,NB*NB*NB*NOB,0.0);
    std::fill_n(I_milj_BB ,NB*NB*NOB*NOB,0.0);
    std::fill_n(I_mibj_BBB,NB*NOB*NOB*NVB,0.0);

    for(auto sg = 0; sg < NB; sg++)
    for(auto lm = 0; lm < NB; lm++)
    for(auto i = 0; i < NOB; i++)
    for(auto nu = 0; nu < NB; nu++)  
    for(auto mu = 0; mu < NB; mu++) {
      I_mils_B[mu + NB*i + NB*NOB*lm + NB*NB*NOB*sg] +=
        (*this->wfn_->moB())(nu,i) *
        (*this->wfn_->aointegrals()->aoERI_)(mu,nu,lm,sg);
    }

    for(auto j = 0; j < NOB; j++)
    for(auto sg = 0; sg < NB; sg++)
    for(auto lm = 0; lm < NB; lm++)
    for(auto i = 0; i < NOB; i++)
    for(auto mu = 0; mu < NB; mu++) {
      I_milj_BB[mu + NB*i + NB*NOB*lm + NB*NB*NOB*j] +=
        (*this->wfn_->moB())(sg,j) *
        I_mils_B[mu + NB*i + NB*NOB*lm + NB*NB*NOB*sg];
    }

    for(auto j = 0; j < NOB; j++)
    for(auto b = 0; b < NVB; b++)
    for(auto lm = 0; lm < NB; lm++)
    for(auto i = 0; i < NOB; i++)
    for(auto mu = 0; mu < NB; mu++) {
      I_mibj_BBB[mu + NB*i + NB*NOB*b + NB*NVB*NOB*j] +=
        (*this->wfn_->moB())(lm,b+NOB) *
        I_milj_BB[mu + NB*i + NB*NOB*lm + NB*NB*NOB*j];
    }

    for(auto j = 0; j < NOB; j++)
    for(auto b = 0; b < NVB; b++)
    for(auto i = 0; i < NOB; i++)
    for(auto a = 0; a < NVB; a++)
    for(auto mu = 0; mu < NB; mu++) {
      VOVOBBBB_[a + NVB*i + NVB*NOB*b + NVB*NVB*NOB*j] +=
        (*this->wfn_->moB())(mu,a+NOB) *
        I_mibj_BBB[mu + NB*i + NB*NOB*b + NB*NVB*NOB*j];
    }


    // Deallocate the unused intermediates
    this->memManager_->free(I_mibj_BBB,NB*NOB*NOB*NVB);
    this->memManager_->free(I_milj_BB ,NB*NB*NOB*NOB);
    this->memManager_->free(I_mils_B  ,NB*NB*NB*NOB);

  }

  // Realize that spin blocks are equivalent for RHF
  if(!isOpenShell) {
    VOVOAABB_ = VOVOAAAA_;
    VOVOBBBB_ = VOVOAAAA_;
  }

  this->haveMOVOVO_ = true;


} // formVOVO

template<>
void MOIntegrals<double>::formFullVOVO(){
  this->formVOVO();
  int NOA = this->wfn_->nOA();
  int NVA = this->wfn_->nVA();
  int NOB = this->wfn_->nOB();
  int NVB = this->wfn_->nVB();
  int NO = this->wfn_->nO();
  int NV = this->wfn_->nV();

  this->VOVO_ = this->memManager_->malloc<double>(NV*NO*NV*NO);
  std::fill_n(this->VOVO_,NO*NV*NO*NV,0.0);

/*
  for(auto J = 0; J < NO; J+=2)
  for(auto B = 0; B < NV; B+=2)
  for(auto I = 0; I < NO; I+=2)
  for(auto A = 0; A < NV; A+=2){
    int a = A / 2;
    int i = I / 2;
    int b = B / 2;
    int j = J / 2;
    this->VOVO_[A + I*NV + B*NO*NV + J*NO*NV*NV] = 
      this->VOVOAAAA_[a + i*NVA + b*NOA*NVA + j*NOA*NVA*NVA];

    this->VOVO_[(A+1) + (I+1)*NV + B*NO*NV + J*NO*NV*NV] = 
      this->VOVOAABB_[a + i*NVA + b*NOA*NVA + j*NOA*NVB*NVA];

    this->VOVO_[A + I*NV + (B+1)*NO*NV + (J+1)*NO*NV*NV] = 
      this->VOVOAABB_[a + i*NVA + b*NOA*NVA + j*NOA*NVB*NVA];

    this->VOVO_[(A+1) + (I+1)*NV + (B+1)*NO*NV + (J+1)*NO*NV*NV] = 
      this->VOVOBBBB_[a + i*NVB + b*NOB*NVB + j*NOB*NVB*NVB];
  }
*/

  for(auto j = 0; j < NOA; j++)
  for(auto b = 0; b < NVA; b++)
  for(auto i = 0; i < NOA; i++)
  for(auto a = 0; a < NVA; a++){
    auto A = 2*a;
    auto I = 2*i;
    auto B = 2*b;
    auto J = 2*j;

    this->VOVO_[A + I*NV + B*NO*NV + J*NO*NV*NV] = 
      this->VOVOAAAA_[a + i*NVA + b*NOA*NVA + j*NOA*NVA*NVA];
  }

  for(auto j = 0; j < NOB; j++)
  for(auto b = 0; b < NVB; b++)
  for(auto i = 0; i < NOA; i++)
  for(auto a = 0; a < NVA; a++){
    auto A = 2*a;
    auto I = 2*i;
    auto B = 2*b + 1;
    auto J = 2*j + 1;

    this->VOVO_[A + I*NV + B*NO*NV + J*NO*NV*NV] = 
      this->VOVOAABB_[a + i*NVA + b*NOA*NVA + j*NOA*NVB*NVA];
  }

/*
  for(auto j = 0; j < NOA; j++)
  for(auto b = 0; b < NVA; b++)
  for(auto i = 0; i < NOB; i++)
  for(auto a = 0; a < NVB; a++){
    auto A = 2*a + 1;
    auto I = 2*i + 1;
    auto B = 2*b;
    auto J = 2*j;

    this->VOVO_[A + I*NV + B*NO*NV + J*NO*NV*NV] = 
      this->VOVOAABB_[b + j*NVA + a*NOA*NVA + i*NOA*NVB*NVA];
  }
*/
/*
  for(auto j = 0; j < NOB; j++)
  for(auto b = 0; b < NVB; b++)
  for(auto i = 0; i < NOA; i++)
  for(auto a = 0; a < NVA; a++){
    auto A = 2*a;
    auto I = 2*i;
    auto B = 2*b + 1;
    auto J = 2*j + 1;

    this->VOVO_[B + J*NV + A*NO*NV + I*NO*NV*NV] = 
      this->VOVOAABB_[a + i*NVA + b*NOA*NVA + j*NOA*NVB*NVA];
  }
*/
  for(auto j = 0; j < NOA; j++)
  for(auto b = 0; b < NVA; b++)
  for(auto i = 0; i < NOB; i++)
  for(auto a = 0; a < NVB; a++){
    auto A = 2*a + 1;
    auto I = 2*i + 1;
    auto B = 2*b;
    auto J = 2*j;

    this->VOVO_[A + I*NV + B*NO*NV + J*NO*NV*NV] =
      this->VOVO_[B + J*NV + A*NO*NV + I*NO*NV*NV];
  }

  for(auto j = 0; j < NOB; j++)
  for(auto b = 0; b < NVB; b++)
  for(auto i = 0; i < NOB; i++)
  for(auto a = 0; a < NVB; a++){
    auto A = 2*a + 1;
    auto I = 2*i + 1;
    auto B = 2*b + 1;
    auto J = 2*j + 1;

    this->VOVO_[A + I*NV + B*NO*NV + J*NO*NV*NV] = 
      this->VOVOBBBB_[a + i*NVB + b*NOB*NVB + j*NOB*NVB*NVB];
  }
} // formFullVOVO

};
