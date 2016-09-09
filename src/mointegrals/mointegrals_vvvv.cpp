#include <mointegrals.h>

namespace ChronusQ {
template<>
void MOIntegrals<double>::formVVVV() {
  if(this->haveMOVVVV_) return;

  bool is2C        = this->wfn_->nTCS() == 2;
  bool isOpenShell = !this->wfn_->isClosedShell;
  int NB           = this->wfn_->nBasis();
  int intDim       = this->pscf_.nVAVA() * this->pscf_.nVAVA();
  int NV = this->wfn_->nVA();
  int NO = this->wfn_->nOA();


  VVVVAAAA_ = this->memManager_->malloc<double>(intDim);
  std::fill_n(VVVVAAAA_,NV*NV*NV*NV,0.0);

  double * I_mnld_A   = this->memManager_->malloc<double>(NB*NB*NB*NV);
  double * I_mncd_AA  = this->memManager_->malloc<double>(NB*NB*NV*NV);
  double * I_mbcd_AAA = this->memManager_->malloc<double>(NB*NV*NV*NV);

  std::fill_n(I_mnld_A  ,NB*NB*NB*NV,0.0);
  std::fill_n(I_mncd_AA ,NB*NB*NV*NV,0.0);
  std::fill_n(I_mbcd_AAA,NB*NV*NV*NV,0.0);

  for(auto d = 0; d < NV; d++)
  for(auto sg = 0; sg < NB; sg++)
  for(auto lm = 0; lm < NB; lm++)
  for(auto nu = 0; nu < NB; nu++)  
  for(auto mu = 0; mu < NB; mu++) {
    I_mnld_A[mu + NB*nu + NB*NB*lm + NB*NB*NB*d] +=
      (*this->wfn_->moA())(sg,d+NO) *
      (*this->wfn_->aointegrals()->aoERI_)(mu,nu,lm,sg);
  }

  for(auto d = 0; d < NV; d++)
  for(auto c = 0; c < NV; c++)
  for(auto lm = 0; lm < NB; lm++)
  for(auto nu = 0; nu < NB; nu++)  
  for(auto mu = 0; mu < NB; mu++) {
    I_mncd_AA[mu + NB*nu + NB*NB*c + NB*NB*NV*d] +=
      (*this->wfn_->moA())(lm,c+NO) *
      I_mnld_A[mu + NB*nu + NB*NB*lm + NB*NB*NB*d];
  }

  for(auto d = 0; d < NV; d++)
  for(auto c = 0; c < NV; c++)
  for(auto b = 0; b < NV; b++)
  for(auto nu = 0; nu < NB; nu++)  
  for(auto mu = 0; mu < NB; mu++) {
    I_mbcd_AAA[mu + NB*b + NB*NV*c + NB*NV*NV*d] +=
      (*this->wfn_->moA())(nu,b+NO) *
      I_mncd_AA[mu + NB*nu + NB*NB*c + NB*NB*NV*d];
  }

  for(auto d = 0; d < NV; d++)
  for(auto c = 0; c < NV; c++)
  for(auto b = 0; b < NV; b++)
  for(auto a = 0; a < NV; a++)
  for(auto mu = 0; mu < NB; mu++) {
    VVVVAAAA_[a + NV*b + NV*NV*c + NV*NV*NV*d] +=
      (*this->wfn_->moA())(mu,a+NO) *
      I_mbcd_AAA[mu + NB*b + NB*NV*c + NB*NV*NV*d];
  }

  this->haveMOVVVV_ = true;

  this->memManager_->free(I_mnld_A  ,NB*NB*NB*NV);
  this->memManager_->free(I_mncd_AA ,NB*NB*NV*NV);
  this->memManager_->free(I_mbcd_AAA,NB*NV*NV*NV);
} // formVVVV

template<>
void MOIntegrals<double>::formFullVVVV(){
  this->formVVVV();
  int NOA = this->wfn_->nOA();
  int NVA = this->wfn_->nVA();
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
      this->VVVVAAAA_[a + b*NVA + c*NVA*NVA + d*NVA*NVA*NVA];

    this->VVVV_[A + B*NV + (C+1)*NV*NV + (D+1)*NV*NV*NV] = 
      this->VVVVAAAA_[a + b*NVA + c*NVA*NVA + d*NVA*NVA*NVA];

    this->VVVV_[(A+1) + (B+1)*NV + (C+1)*NV*NV + (D+1)*NV*NV*NV] = 
      this->VVVVAAAA_[a + b*NVA + c*NVA*NVA + d*NVA*NVA*NVA];
  }
} // formFullVVVV

};
