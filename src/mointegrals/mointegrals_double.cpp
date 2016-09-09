#include <mointegrals.h>

namespace ChronusQ {
template<>
void MOIntegrals<double>::getLocalMOs() {
  if(this->haveLocMO_) return;
  bool is2C        = this->wfn_->nTCS() == 2;
  bool isOpenShell = !this->wfn_->isClosedShell;
  int NB           = this->wfn_->nBasis();
//int NBT          = this->wfn_->nTCS() * NB;


/*
  if(is2C) {
    this->locMOAOcc_ = 
      std::unique_ptr<RealTensor2d>(new RealTensor2d(NB,this->wfn_->nO()));
    this->locMOBOcc_ = 
      std::unique_ptr<RealTensor2d>(new RealTensor2d(NB,this->wfn_->nO()));
    this->locMOAVir_ = 
      std::unique_ptr<RealTensor2d>(new RealTensor2d(NB,this->wfn_->nV()));
    this->locMOBVir_ = 
      std::unique_ptr<RealTensor2d>(new RealTensor2d(NB,this->wfn_->nV()));

    for(auto i  = 0; i  < this->wfn_->nO(); i++ )
    for(auto mu = 0; mu < NB              ; mu++){
      (*this->locMOAOcc_)(mu,i) = (*this->wfn_->moA())(2*mu    ,i);
      (*this->locMOBOcc_)(mu,i) = (*this->wfn_->moA())(2*mu + 1,i);
    }
    for(auto a  = 0; a  < this->wfn_->nV(); a++ )
    for(auto mu = 0; mu < NB              ; mu++){
      int A = this->wfn_->nO() + a;       
      (*this->locMOAVir_)(mu,a) = (*this->wfn_->moA())(2*mu    ,A);
      (*this->locMOBVir_)(mu,a) = (*this->wfn_->moA())(2*mu + 1,A);
    }
  } else {
    this->locMOAOcc_ = 
      std::unique_ptr<RealTensor2d>(new RealTensor2d(NB,this->wfn_->nOA()));
    this->locMOAVir_ = 
      std::unique_ptr<RealTensor2d>(new RealTensor2d(NB,this->wfn_->nVA()));

    std::copy(
      this->wfn_->moA()->data(),
      this->wfn_->moA()->data() + this->wfn_->nOA(),
      &this->locMOAOcc_->storage()[0]);
    std::copy(
      this->wfn_->moA()->data() + this->wfn_->nOA(),
      this->wfn_->moA()->data() + this->wfn_->moA()->size(),
      &this->locMOAVir_->storage()[0]);

    if(isOpenShell) {
      this->locMOBOcc_ = 
        std::unique_ptr<RealTensor2d>(new RealTensor2d(NB,this->wfn_->nOB()));
      this->locMOBVir_ = 
        std::unique_ptr<RealTensor2d>(new RealTensor2d(NB,this->wfn_->nVB()));
      std::copy(
        this->wfn_->moB()->data(),
        this->wfn_->moB()->data() + this->wfn_->nOB(),
        &this->locMOBOcc_->storage()[0]);
      std::copy(
        this->wfn_->moB()->data() + this->wfn_->nOB(),
        this->wfn_->moB()->data() + this->wfn_->moB()->size(),
        &this->locMOBVir_->storage()[0]);
    }
  }
*/

}; // getLocalMOs

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

}

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
}

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
}

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
}

template<>
void MOIntegrals<double>::testMOInts(){
  double EMP2 = 0;
  int NO = this->wfn_->nO();
  int NV = this->wfn_->nV();
  int NOA = this->wfn_->nOA();
  int NVA = this->wfn_->nVA();
  int NOB = this->wfn_->nOB();
  int NVB = this->wfn_->nVB();

  double * EPSA = this->wfn_->epsA()->data();
  double * EPSB = this->wfn_->isClosedShell ? EPSA : this->wfn_->epsB()->data();

  this->formFullVOVO();

/*
  // Spacial orbital MP2 energy
  for(auto j = 0; j < NO; j++)
  for(auto i = 0; i < NO; i++)
  for(auto b = 0; b < NV; b++)
  for(auto a = 0; a < NV; a++){
    EMP2 += VOVOAAAA_[a + NV*i + NV*NO*b + NV*NV*NO*j] *
      ( 2*VOVOAAAA_[a + NV*i + NV*NO*b + NV*NV*NO*j]
          - VOVOAAAA_[b + NV*i + NV*NO*a + NV*NV*NO*j]) /
    ((*this->wfn_->epsA())(i) + (*this->wfn_->epsA())(j)
    -(*this->wfn_->epsA())(a+NO) - (*this->wfn_->epsA())(b+NO));
  }
  cout << "EMP2 = " << EMP2 << endl;
*/

  EMP2 = 0;
  double eI, eJ, eA, eB;
  for(int j = 0; j < NO; j++)
  for(int i = 0; i < NO; i++)
  for(int b = 0; b < NV; b++)
  for(int a = 0; a < NV; a++){
    
/*
    double eI = (*this->wfn_->epsA())(i/2);
    double eJ = (*this->wfn_->epsA())(j/2);
    double eA = (*this->wfn_->epsA())(a/2 + NOA);
    double eB = (*this->wfn_->epsA())(b/2 + NOA);
*/
    eI = (i % 2 == 0) ? EPSA[i/2]       : EPSB[i/2];
    eJ = (j % 2 == 0) ? EPSA[j/2]       : EPSB[j/2];
    eA = (a % 2 == 0) ? EPSA[a/2 + NOA] : EPSB[a/2 + NOB];
    eB = (b % 2 == 0) ? EPSA[b/2 + NOA] : EPSB[b/2 + NOB];


    // < IJ || AB > = (AI | BJ) - (BI | AJ)
    double DiracIJAB = 
      VOVO_[a + i*NV + b*NO*NV + j*NO*NV*NV] -
      VOVO_[b + i*NV + a*NO*NV + j*NO*NV*NV];

    // dIJAB = e(I) + e(J) - e(A) - e(B)
    double deltaIJAB = eI + eJ - eA - eB;

    // E(2) += |< IJ || AB >|^2 / dIJAB
    EMP2 += DiracIJAB * DiracIJAB / deltaIJAB;
  }
  
  // E(2) = E(2) / 4
  EMP2 *= 0.25;
  cout << "EMP2 = " << EMP2 << endl;


/*
  this->formFullVVOO();
  this->formFullVVVV();
  this->formFullOOOO();

  double EMP3_1 = 0;
  double EMP3_2 = 0;
  double EMP3_3 = 0;

  for(auto i = 0; i < 2*NO; i++)
  for(auto j = 0; j < 2*NO; j++)
  for(auto k = 0; k < 2*NO; k++)
  for(auto l = 0; l < 2*NO; l++)
  for(auto a = 0; a < 2*NV; a++)
  for(auto b = 0; b < 2*NV; b++){

    double eI = (*this->wfn_->epsA())(i/2);
    double eJ = (*this->wfn_->epsA())(j/2);
    double eK = (*this->wfn_->epsA())(k/2);
    double eL = (*this->wfn_->epsA())(l/2);
    double eA = (*this->wfn_->epsA())(a/2 + NOA);
    double eB = (*this->wfn_->epsA())(b/2 + NOA);

    // < IJ || AB > = (AI | BJ) - (BI | AJ)
    double DiracIJAB = 
      VOVO_[a + i*2*NV + b*4*NO*NV + j*8*NV*NV*NO] - 
      VOVO_[b + i*2*NV + a*4*NO*NV + j*8*NV*NV*NO]; 

    // < IJ || KL > = (IK | JL) - (IL | JK)
    double DiracIJKL = 
      OOOO_[i + k*2*NO + j*4*NO*NO + l*8*NO*NO*NO] - 
      OOOO_[i + l*2*NO + j*4*NO*NO + k*8*NO*NO*NO]; 

    // < KL || AB > = (AK | BL) - (BK | AL)
    double DiracKLAB = 
      VOVO_[a + k*2*NV + b*4*NO*NV + l*8*NV*NV*NO] - 
      VOVO_[b + k*2*NV + a*4*NO*NV + l*8*NV*NV*NO]; 


    // dIJAB = e(I) + e(J) - e(A) - e(B)
    double deltaIJAB = eI + eJ - eA - eB;

    // dKLAB = e(K) + e(L) - e(A) - e(B)
    double deltaKLAB= eK + eL - eA - eB;

    // E(3,1) += < IJ || AB > * < IJ || KL > * < KL || AB > / dIJAB / dKLAB 
    EMP3_1 += DiracIJAB * DiracIJKL * DiracKLAB / deltaIJAB / deltaKLAB;
  }

  for(auto i = 0; i < 2*NO; i++)
  for(auto j = 0; j < 2*NO; j++)
  for(auto a = 0; a < 2*NV; a++)
  for(auto b = 0; b < 2*NV; b++) 
  for(auto c = 0; c < 2*NV; c++)
  for(auto d = 0; d < 2*NV; d++){

    double eI = (*this->wfn_->epsA())(i/2);
    double eJ = (*this->wfn_->epsA())(j/2);
    double eA = (*this->wfn_->epsA())(a/2 + NOA);
    double eB = (*this->wfn_->epsA())(b/2 + NOA);
    double eC = (*this->wfn_->epsA())(c/2 + NOA);
    double eD = (*this->wfn_->epsA())(d/2 + NOA);

    // < IJ || AB > = (AI | BJ) - (BI | AJ)
    double DiracIJAB = 
      VOVO_[a + i*2*NV + b*4*NO*NV + j*8*NV*NV*NO] - 
      VOVO_[b + i*2*NV + a*4*NO*NV + j*8*NV*NV*NO]; 

    // < AB || CD > = (AC | BD) - (AD | BC)
    double DiracABCD = 
      VVVV_[a + c*2*NV + b*4*NV*NV + d*8*NV*NV*NV] - 
      VVVV_[a + d*2*NV + b*4*NV*NV + c*8*NV*NV*NV]; 

    // < IJ || CD > = (CI | DJ) - (DI | CJ)
    double DiracIJCD = 
      VOVO_[c + i*2*NV + d*4*NO*NV + j*8*NV*NV*NO] - 
      VOVO_[d + i*2*NV + c*4*NO*NV + j*8*NV*NV*NO]; 


    // dIJAB = e(I) + e(J) - e(A) - e(B)
    double deltaIJAB = eI + eJ - eA - eB;

    // dIJCD = e(I) + e(J) - e(C) - e(D)
    double deltaIJCD = eI + eJ - eC - eD; 

    // E(3,2) += < IJ || AB > * < AB || CD > * < IJ || CD > / dIJAB / dIJCD 
    EMP3_2 += DiracIJAB * DiracABCD * DiracIJCD / deltaIJAB / deltaIJCD;
    
  }

  for(auto i = 0; i < 2*NO; i++)
  for(auto j = 0; j < 2*NO; j++)
  for(auto k = 0; k < 2*NO; k++)
  for(auto a = 0; a < 2*NV; a++)
  for(auto b = 0; b < 2*NV; b++) 
  for(auto c = 0; c < 2*NV; c++){

    double eI = (*this->wfn_->epsA())(i/2);
    double eJ = (*this->wfn_->epsA())(j/2);
    double eK = (*this->wfn_->epsA())(k/2);
    double eA = (*this->wfn_->epsA())(a/2 + NOA);
    double eB = (*this->wfn_->epsA())(b/2 + NOA);
    double eC = (*this->wfn_->epsA())(c/2 + NOA);

    // < IJ || AB > = (AI | BJ) - (BI | AJ)
    double DiracIJAB = 
      VOVO_[a + i*2*NV + b*4*NO*NV + j*8*NV*NV*NO] - 
      VOVO_[b + i*2*NV + a*4*NO*NV + j*8*NV*NV*NO]; 

    // < BK || CJ > = (BC | KJ) - (BJ | CK)
    double DiracBKCJ = 
      VVOO_[b + c*2*NV + k*4*NV*NV + j*8*NV*NV*NO] -
      VOVO_[b + j*2*NV + c*4*NO*NV + k*8*NO*NV*NV];

    // < IK || AC > = (AI | CK) - (CI | AK)
    double DiracIKAC = 
      VOVO_[a + i*2*NV + c*4*NO*NV + k*8*NV*NV*NO] - 
      VOVO_[c + i*2*NV + a*4*NO*NV + k*8*NV*NV*NO]; 


    // dIJAB = e(I) + e(J) - e(A) - e(B)
    double deltaIJAB = eI + eJ - eA - eB;

    // dIKAC = e(I) + e(K) - e(A) - e(C)
    double deltaIKAC = eI + eK - eA - eC; 

    // E(3,3) -= < IJ || AB > * < BK || CJ > * < IK || AC > / dIJAB / dIKAC
    EMP3_3 -= DiracIJAB * DiracBKCJ * DiracIKAC / deltaIJAB / deltaIKAC;
  }

  // E(3) = (1/8) * ( E(3,1) + E(3,2) ) + E(3,3)
  double EMP3 = (0.125)*(EMP3_1 + EMP3_2) + EMP3_3;

  cout << "EMP3 = " << EMP3 << endl;
*/
};
};
