#include <mointegrals.h>

namespace ChronusQ {
template<>
void MOIntegrals<double>::formVVOO() {
  if(this->haveMOVVOO_) return;

  bool is2C        = this->wfn_->nTCS() == 2;
  bool isOpenShell = !this->wfn_->isClosedShell;
  int NB           = this->wfn_->nBasis();
  int intDim       = this->pscf_.nOAVA() * this->pscf_.nOAVA();
  int NO = this->wfn_->nOA();
  int NV = this->wfn_->nVA();


  VVOOAAAA_ = this->memManager_->malloc<double>(intDim);
  std::fill_n(VVOOAAAA_,NV*NO*NO*NV,0.0);

  double * I_mnis_A   = this->memManager_->malloc<double>(NB*NB*NB*NO);
  double * I_mnij_AA  = this->memManager_->malloc<double>(NB*NB*NO*NO);
  double * I_mbij_AAA = this->memManager_->malloc<double>(NB*NO*NO*NV);

  std::fill_n(I_mnis_A  ,NB*NB*NB*NO,0.0);
  std::fill_n(I_mnij_AA ,NB*NB*NO*NO,0.0);
  std::fill_n(I_mbij_AAA,NB*NO*NO*NV,0.0);

  for(auto sg = 0; sg < NB; sg++)
  for(auto i = 0; i < NO; i++)
  for(auto lm = 0; lm < NB; lm++)
  for(auto nu = 0; nu < NB; nu++)  
  for(auto mu = 0; mu < NB; mu++) {
    I_mnis_A[mu + NB*nu + NB*NB*i + NB*NB*NO*sg] +=
      (*this->wfn_->moA())(lm,i) *
      (*this->wfn_->aointegrals()->aoERI_)(mu,nu,lm,sg);
  }

  for(auto j = 0; j < NO; j++)
  for(auto sg = 0; sg < NB; sg++)
  for(auto i = 0; i < NO; i++)
  for(auto nu = 0; nu < NB; nu++)  
  for(auto mu = 0; mu < NB; mu++) {
    I_mnij_AA[mu + NB*nu + NB*NB*i + NB*NB*NO*j] +=
      (*this->wfn_->moA())(sg,j) *
      I_mnis_A[mu + NB*nu + NB*NB*i + NB*NB*NO*sg];
  }

  for(auto j = 0; j < NO; j++)
  for(auto i = 0; i < NO; i++)
  for(auto b = 0; b < NV; b++)
  for(auto nu = 0; nu < NB; nu++)  
  for(auto mu = 0; mu < NB; mu++) {
    I_mbij_AAA[mu + NB*b + NV*NB*i + NB*NV*NO*j] +=
      (*this->wfn_->moA())(nu,b+NO) *
      I_mnij_AA[mu + NB*nu + NB*NB*i + NB*NB*NO*j];
  }

  for(auto j = 0; j < NO; j++)
  for(auto i = 0; i < NO; i++)
  for(auto b = 0; b < NV; b++)
  for(auto a = 0; a < NV; a++)
  for(auto mu = 0; mu < NB; mu++) {
    VVOOAAAA_[a + NV*b + NV*NV*i + NV*NV*NO*j] +=
      (*this->wfn_->moA())(mu,a+NO) *
      I_mbij_AAA[mu + NB*b + NV*NB*i + NB*NV*NO*j];
  }

  this->haveMOVVOO_ = true;

  this->memManager_->free(I_mnis_A  ,NB*NB*NB*NO);
  this->memManager_->free(I_mnij_AA ,NB*NB*NO*NO);
  this->memManager_->free(I_mbij_AAA,NB*NO*NO*NV);
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
