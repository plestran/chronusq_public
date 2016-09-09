#include <mointegrals.h>

namespace ChronusQ {
template<>
void MOIntegrals<double>::formOOOO() {
  if(this->haveMOOOOO_) return;

  bool is2C        = this->wfn_->nTCS() == 2;
  bool isOpenShell = !this->wfn_->isClosedShell;
  int NB           = this->wfn_->nBasis();
  int intDim       = this->pscf_.nOAOA() * this->pscf_.nOAOA();
  int NV = this->wfn_->nVA();
  int NO = this->wfn_->nOA();


  OOOOAAAA_ = this->memManager_->malloc<double>(intDim);
  std::fill_n(OOOOAAAA_,NO*NO*NO*NO,0.0);

  double * I_1_A   = this->memManager_->malloc<double>(NB*NB*NB*NO);
  double * I_2_AA  = this->memManager_->malloc<double>(NB*NB*NO*NO);
  double * I_3_AAA = this->memManager_->malloc<double>(NB*NO*NO*NO);

  std::fill_n(I_1_A  ,NB*NB*NB*NO,0.0);
  std::fill_n(I_2_AA ,NB*NB*NO*NO,0.0);
  std::fill_n(I_3_AAA,NB*NO*NO*NO,0.0);

  for(auto l = 0; l < NO; l++)
  for(auto sg = 0; sg < NB; sg++)
  for(auto lm = 0; lm < NB; lm++)
  for(auto nu = 0; nu < NB; nu++)  
  for(auto mu = 0; mu < NB; mu++) {
    I_1_A[mu + NB*nu + NB*NB*lm + NB*NB*NB*l] +=
      (*this->wfn_->moA())(sg,l) *
      (*this->wfn_->aointegrals()->aoERI_)(mu,nu,lm,sg);
  }

  for(auto l = 0; l < NO; l++)
  for(auto k = 0; k < NO; k++)
  for(auto lm = 0; lm < NB; lm++)
  for(auto nu = 0; nu < NB; nu++)  
  for(auto mu = 0; mu < NB; mu++) {
    I_2_AA[mu + NB*nu + NB*NB*k + NB*NB*NO*l] +=
      (*this->wfn_->moA())(lm,k) *
      I_1_A[mu + NB*nu + NB*NB*lm + NB*NB*NB*l];
  }

  for(auto l = 0; l < NO; l++)
  for(auto k = 0; k < NO; k++)
  for(auto j = 0; j < NO; j++)
  for(auto nu = 0; nu < NB; nu++)  
  for(auto mu = 0; mu < NB; mu++) {
    I_3_AAA[mu + NB*j + NB*NO*k + NB*NO*NO*l] +=
      (*this->wfn_->moA())(nu,j) *
      I_2_AA[mu + NB*nu + NB*NB*k + NB*NB*NO*l];
  }

  for(auto l = 0; l < NO; l++)
  for(auto k = 0; k < NO; k++)
  for(auto j = 0; j < NO; j++)
  for(auto i = 0; i < NO; i++)
  for(auto mu = 0; mu < NB; mu++) {
    OOOOAAAA_[i + NO*j + NO*NO*k + NO*NO*NO*l] +=
      (*this->wfn_->moA())(mu,i) *
      I_3_AAA[mu + NB*j + NB*NO*k + NB*NO*NO*l];
  }

  this->haveMOOOOO_ = true;

  this->memManager_->free(I_1_A  ,NB*NB*NB*NO);
  this->memManager_->free(I_2_AA ,NB*NB*NO*NO);
  this->memManager_->free(I_3_AAA,NB*NO*NO*NO);
} // formOOOO

template<>
void MOIntegrals<double>::formFullOOOO(){
  this->formOOOO();
  int NOA = this->wfn_->nOA();
  int NVA = this->wfn_->nVA();
  int NO = this->wfn_->nO();
  int NV = this->wfn_->nV();

  this->OOOO_ = this->memManager_->malloc<double>(NO*NO*NO*NO);
  std::fill_n(this->OOOO_,NO*NO*NO*NO,0.0);

  for(auto L = 0; L < NO; L+=2)
  for(auto K = 0; K < NO; K+=2)
  for(auto J = 0; J < NO; J+=2)
  for(auto I = 0; I < NO; I+=2){
    int i = I / 2;
    int j = J / 2;
    int k = K / 2;
    int l = L / 2;

    this->OOOO_[I + J*NO + K*NO*NO + L*NO*NO*NO] = 
      this->OOOOAAAA_[i + j*NOA + k*NOA*NOA + l*NOA*NOA*NOA];

    this->OOOO_[(I+1) + (J+1)*NO + K*NO*NO + L*NO*NO*NO] = 
      this->OOOOAAAA_[i + j*NOA + k*NOA*NOA + l*NOA*NOA*NOA];

    this->OOOO_[I + J*NO + (K+1)*NO*NO + (L+1)*NO*NO*NO] = 
      this->OOOOAAAA_[i + j*NOA + k*NOA*NOA + l*NOA*NOA*NOA];

    this->OOOO_[(I+1) + (J+1)*NO + (K+1)*NO*NO + (L+1)*NO*NO*NO] = 
      this->OOOOAAAA_[i + j*NOA + k*NOA*NOA + l*NOA*NOA*NOA];
  }
} // formFullOOOO

};
