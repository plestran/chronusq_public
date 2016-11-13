/* 
 *  This file is part of the Chronus Quantum (ChronusQ) software package
 *  
 *  Copyright (C) 2014-2016 Li Research Group (University of Washington)
 *  
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *  
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *  
 *  You should have received a copy of the GNU General Public License along
 *  with this program; if not, write to the Free Software Foundation, Inc.,
 *  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *  
 *  Contact the Developers:
 *    E-Mail: xsli@uw.edu
 *  
 */
template <typename T>
void MOIntegrals<T>::formVOVO() {
  if(this->haveMOVOVO_) return;

  bool is2C        = this->wfn_->nTCS() == 2;
  bool isOpenShell = !this->wfn_->isClosedShell;
  int NB           = this->wfn_->nBasis();
  int NOA = this->wfn_->nOA();
  int NVA = this->wfn_->nVA();
  int NOB = this->wfn_->nOB();
  int NVB = this->wfn_->nVB();

  if(is2C) {
    form2CVOVO();
    return;
  }

/*
  if(isOpenShell) {
    *this->wfn_->moB() = *this->wfn_->moA();
  //this->wfn_->moB()->col(1) *= -1;
    prettyPrint(cout,*this->wfn_->moA(),"MOA");
    prettyPrint(cout,*this->wfn_->moB(),"MOB");
  }
*/

  /// Evaluate VOVO (AAAA) MO Integrals

  // Allocate and zero out VOVO (AAAA) storage
  VOVOAAAA_ = this->memManager_->template malloc<T>(NOA*NVA*NOA*NVA);
  std::fill_n(VOVOAAAA_,NVA*NOA*NOA*NVA,0.0);

  // Allocate and zero out Itermediates
  T * I_mils_A   = this->memManager_->template malloc<T>(NB*NB*NB*NOA);
  T * I_milj_AA  = this->memManager_->template malloc<T>(NB*NB*NOA*NOA);
  T * I_mibj_AAA = this->memManager_->template malloc<T>(NB*NOA*NOA*NVA);

  std::fill_n(I_mils_A  ,NB*NB*NB*NOA,0.0);
  std::fill_n(I_milj_AA ,NB*NB*NOA*NOA,0.0);
  std::fill_n(I_mibj_AAA,NB*NOA*NOA*NVA,0.0);

  // First Quarter (AAAA) transformation (mn|ls) -> (mi|ls) (A)
  rank4w2Contract(2,this->wfn_->moA()->data(),NB,
    &this->wfn_->aointegrals()->aoERI_->storage()[0],NB,NB,NB,NB,
    I_mils_A,NOA);

  // Second Quarter (AAAA) transformation (mi|ls) (A) -> (mi|lj) (AA)
  rank4w2Contract(4,this->wfn_->moA()->data(),NB,
    I_mils_A,NB,NOA,NB,NB,
    I_milj_AA,NOA);

  // Conjugate for contraction if complex
  if(std::is_same<dcomplex,T>::value) {
    (*this->wfn_->moA()) = this->wfn_->moA()->conjugate();
  }

  // Third Quarter (AAAA) transformation (mi|lj) (AA) -> (mi|bj) (AAA)
  rank4w2Contract(3,this->wfn_->moA()->data() + NB*NOA,NB,
    I_milj_AA,NB,NOA,NB,NOA,
    I_mibj_AAA,NVA);

  // Fourth Quarter (AAAA) transformation (mi|bj) (AAA) -> (ai|bj) (AAAA)
  rank4w2Contract(1,this->wfn_->moA()->data() + NB*NOA,NB,
    I_mibj_AAA,NB,NOA,NVA,NOA,
    VOVOAAAA_,NVA);

  // Conjugate Back for contraction if complex
  if(std::is_same<dcomplex,T>::value) {
    (*this->wfn_->moA()) = this->wfn_->moA()->conjugate();
  }


  // Deallocate the unused intermediates
  this->memManager_->free(I_mibj_AAA,NB*NOA*NOA*NVA);
  this->memManager_->free(I_milj_AA ,NB*NB*NOA*NOA);

  if(isOpenShell) {
    /// Evaluate VOVO (AABB) MO Integrals

    // Allocate and zero out VOVO (AABB) storage
    VOVOAABB_ = this->memManager_->template malloc<T>(NOA*NVA*NOB*NVB);
    std::fill_n(VOVOAABB_,NOA*NVA*NOB*NVB,0.0);

    // Allocate and zero out Itermediates
    T * I_milj_AB  = this->memManager_->template malloc<T>(NB*NB*NOA*NOB);
    T * I_mibj_ABB = this->memManager_->template malloc<T>(NB*NOA*NOB*NVB);

    std::fill_n(I_milj_AB ,NB*NB*NOA*NOB,0.0);
    std::fill_n(I_mibj_ABB,NB*NOA*NOB*NVB,0.0);

    // Second Quarter (AABB) transformation (mi|ls) (A) -> (mi|lj) (AB)
    rank4w2Contract(4,this->wfn_->moB()->data(),NB,
      I_mils_A,NB,NOA,NB,NB,
      I_milj_AB,NOB);
 
    // Conjugate for contraction if complex
    if(std::is_same<dcomplex,T>::value) {
      (*this->wfn_->moA()) = this->wfn_->moA()->conjugate();
      (*this->wfn_->moB()) = this->wfn_->moB()->conjugate();
    }

    // Third Quarter (AABB) transformation (mi|lj) (AB) -> (mi|bj) (ABB)
    rank4w2Contract(3,this->wfn_->moB()->data() + NB*NOB,NB,
      I_milj_AB,NB,NOA,NB,NOB,
      I_mibj_ABB,NVB);

    // Fourth Quarter (AABB) transformation (mi|bj) (ABB) -> (ai|bj) (AABB)
    rank4w2Contract(1,this->wfn_->moA()->data() + NB*NOA,NB,
      I_mibj_ABB,NB,NOA,NVB,NOB,
      VOVOAABB_,NVA);

    // Conjugate back for contraction if complex
    if(std::is_same<dcomplex,T>::value) {
      (*this->wfn_->moA()) = this->wfn_->moA()->conjugate();
      (*this->wfn_->moB()) = this->wfn_->moB()->conjugate();
    }

    // Deallocate the unused intermediates
    this->memManager_->free(I_milj_AB ,NB*NB*NOA*NOB);
    this->memManager_->free(I_mibj_ABB,NB*NOA*NOB*NVB);
  }

  // Deallocate the unused intermediates
  this->memManager_->free(I_mils_A  ,NB*NB*NB*NOA);

  if(isOpenShell) {
    /// Evaluate VOVO (BBBB) MO Integrals

    // Allocate and zero out VOVO (BBBB) storage
    VOVOBBBB_ = this->memManager_->template malloc<T>(NOB*NVB*NOB*NVB);
    std::fill_n(VOVOBBBB_,NVB*NOB*NOB*NVB,0.0);

    T * I_mils_B   = this->memManager_->template malloc<T>(NB*NB*NB*NOB);
    T * I_milj_BB  = this->memManager_->template malloc<T>(NB*NB*NOB*NOB);
    T * I_mibj_BBB = this->memManager_->template malloc<T>(NB*NOB*NOB*NVB);

    // Allocate and zero out Itermediates
    std::fill_n(I_mils_B  ,NB*NB*NB*NOB,0.0);
    std::fill_n(I_milj_BB ,NB*NB*NOB*NOB,0.0);
    std::fill_n(I_mibj_BBB,NB*NOB*NOB*NVB,0.0);


    // First Quarter (BBBB) transformation (mn|ls) -> (mi|ls) (B)
    rank4w2Contract(2,this->wfn_->moB()->data(),NB,
      &this->wfn_->aointegrals()->aoERI_->storage()[0],NB,NB,NB,NB,
      I_mils_B,NOB);

    // Second Quarter (BBBB) transformation (mi|ls) (B) -> (mi|lj) (BB)
    rank4w2Contract(4,this->wfn_->moB()->data(),NB,
      I_mils_B,NB,NOB,NB,NB,
      I_milj_BB,NOB);

    // Conjugate for contraction if complex
    if(std::is_same<dcomplex,T>::value) {
      (*this->wfn_->moB()) = this->wfn_->moB()->conjugate();
    }

    // Third Quarter (BBBB) transformation (mi|lj) (BB) -> (mi|bj) (BBB)
    rank4w2Contract(3,this->wfn_->moB()->data() + NB*NOB,NB,
      I_milj_BB,NB,NOB,NB,NOB,
      I_mibj_BBB,NVB);

    // Fourth Quarter (BBBB) transformation (mi|bj) (BBB) -> (ai|bj) (BBBB)
    rank4w2Contract(1,this->wfn_->moB()->data() + NB*NOB,NB,
      I_mibj_BBB,NB,NOB,NVB,NOB,
      VOVOBBBB_,NVB);

    // Conjugate back for contraction if complex
    if(std::is_same<dcomplex,T>::value) {
      (*this->wfn_->moB()) = this->wfn_->moB()->conjugate();
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

template <typename T> 
void MOIntegrals<T>::form2CVOVO() {

  int NB = this->wfn_->nBasis();
  int NO = this->wfn_->nO();
  int NV = this->wfn_->nV();

  VOVO_ = this->memManager_->template malloc<T>(NO*NV*NO*NV);
  std::fill_n(VOVO_,NO*NV*NO*NV,0.);

  T* tmp1 = this->memManager_->template malloc<T>(NO*NV*NB*NB);
  std::fill_n(tmp1,NO*NV*NB*NB,0.0);

/*
  // N^8 Algorithm
  for(auto j = 0; j < NO; j++)
  for(auto b = 0; b < NV; b++)
  for(auto i = 0; i < NO; i++)
  for(auto a = 0; a < NV; a++) {

  for(auto sg = 0; sg < 2*NB; sg+=2)
  for(auto lm = 0; lm < 2*NB; lm+=2)
  for(auto nu = 0; nu < 2*NB; nu+=2)
  for(auto mu = 0; mu < 2*NB; mu+=2){

    VOVO_[a + i*NV + b*NO*NV + j*NV*NO*NV] +=
      (
      std::conj((*wfn_->moA())(mu,a+NO)) * (*wfn_->moA())(nu,i) *
      std::conj((*wfn_->moA())(lm,b+NO)) * (*wfn_->moA())(sg,j) +

      std::conj((*wfn_->moA())(mu+1,a+NO)) * (*wfn_->moA())(nu+1,i) *
      std::conj((*wfn_->moA())(lm,b+NO))   * (*wfn_->moA())(sg,j) +

      std::conj((*wfn_->moA())(mu,a+NO))   * (*wfn_->moA())(nu,i) *
      std::conj((*wfn_->moA())(lm+1,b+NO)) * (*wfn_->moA())(sg+1,j) +

      std::conj((*wfn_->moA())(mu+1,a+NO)) * (*wfn_->moA())(nu+1,i) *
      std::conj((*wfn_->moA())(lm+1,b+NO)) * (*wfn_->moA())(sg+1,j) 
      ) * (*wfn_->aointegrals()->aoERI_)(mu/2,nu/2,lm/2,sg/2);
  }
  }
*/



  int nDo = 200;
 
  TMap scr(this->memManager_->template malloc<T>(nDo*NB*NB),NB*NB,nDo);
  TMap scr2(this->memManager_->template malloc<T>(4*NB*NB),2*NB,2*NB);
  RealMap aoERI(&this->wfn_->aointegrals()->aoERI_->storage()[0],NB*NB,NB*NB);
  Eigen::Map<TMatrix,0,Eigen::Stride<Dynamic,Dynamic> > 
    SAA(scr2.data(),NB,NB,Eigen::Stride<Dynamic,Dynamic>(4*NB,2));
  Eigen::Map<TMatrix,0,Eigen::Stride<Dynamic,Dynamic> > 
    SBB(scr2.data()+2*NB+1,NB,NB,Eigen::Stride<Dynamic,Dynamic>(4*NB,2));

  std::vector<TMap> scrMaps;
  for(auto iDo = 0; iDo < nDo; iDo++) 
    scrMaps.emplace_back(scr.data()+iDo*NB*NB,NB,NB);

  this->wfn_->fileio()->out << "Begining First Half Transformation AIBJ"
    << endl;
  for(auto bj = 0; bj < NV*NO; bj+=nDo) {
    auto NDo = std::min(nDo, NO*NV-bj);
    for(auto iDo = 0, BJ = bj; iDo < NDo; iDo++, BJ++){ 
      int b = BJ % NV;
      int j = BJ / NV;
      scr2.noalias() = 
        wfn_->moA()->col(b+NO).conjugate() * wfn_->moA()->col(j).transpose();
      scrMaps[iDo].noalias() = SAA + SBB;
    }

    TMap TMPMAP(tmp1+bj*NB*NB,NB*NB,NDo);
    TMPMAP.noalias() = aoERI * scr.block(0,0,NB*NB,NDo);
  }
 
  this->wfn_->fileio()->out << "Begining Second Half Transformation AIBJ"
    << endl;
  TMap TMPMAP(tmp1,NB*NB,NO*NV);
  TMap VOVOMAP(VOVO_,NO*NV,NO*NV);
  for(auto ai = 0; ai < NV*NO; ai+=nDo) {
    auto NDo = std::min(nDo, NO*NV - ai);
    for(auto iDo = 0, AI = ai; iDo < NDo; iDo++, AI++){ 
      int a = AI % NV;
      int i = AI / NV;
      scr2.noalias() = 
        wfn_->moA()->col(a+NO).conjugate() * wfn_->moA()->col(i).transpose();
      scrMaps[iDo].noalias() = SAA + SBB;
    }

    VOVOMAP.block(ai,0,NDo,NO*NV).noalias() = 
      scr.block(0,0,NB*NB,NDo).transpose() * TMPMAP;
  }


  this->memManager_->free(tmp1,NO*NV*NB*NB);
  this->memManager_->free(scr.data(),nDo*NB*NB);
  this->memManager_->free(scr2.data(),4*NB*NB);

  this->haveMOVOVO_ = true;
}; // form2CVOVO

template <typename T>
void MOIntegrals<T>::formFullVOVO(){
  this->formVOVO();
  if(this->wfn_->nTCS() == 2) return;

  int NOA = this->wfn_->nOA();
  int NVA = this->wfn_->nVA();
  int NOB = this->wfn_->nOB();
  int NVB = this->wfn_->nVB();
  int NO = this->wfn_->nO();
  int NV = this->wfn_->nV();

  this->VOVO_ = this->memManager_->template malloc<T>(NV*NO*NV*NO);
  std::fill_n(this->VOVO_,NO*NV*NO*NV,0.0);

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

  for(auto j = 0; j < NOA; j++)
  for(auto b = 0; b < NVA; b++)
  for(auto i = 0; i < NOB; i++)
  for(auto a = 0; a < NVB; a++){
    auto A = 2*a + 1;
    auto I = 2*i + 1;
    auto B = 2*b;
    auto J = 2*j;

    this->VOVO_[A + I*NV + B*NO*NV + J*NO*NV*NV] =
      this->VOVO_[b + j*NV + a*NO*NV + i*NO*NV*NV];
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
