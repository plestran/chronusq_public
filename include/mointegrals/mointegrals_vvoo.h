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
void MOIntegrals<T>::formVVOO() {
  if(this->haveMOVVOO_) return;

  bool is2C        = this->wfn_->nTCS() == 2;
  bool isOpenShell = !this->wfn_->isClosedShell;
  int NB           = this->wfn_->nBasis();
  int intDim       = this->pscf_.nOAVA() * this->pscf_.nOAVA();
  int NOA = this->wfn_->nOA();
  int NVA = this->wfn_->nVA();
  int NOB = this->wfn_->nOB();
  int NVB = this->wfn_->nVB();

  if(is2C) {
    form2CVVOO();
    return;
  }

  /// Evaluate VVOO (AAAA) MO Integrals

  // Allocate and zero out VVOO (AAAA) storage
  VVOOAAAA_ = this->memManager_->template malloc<T>(NVA*NOA*NOA*NVA);
  std::fill_n(VVOOAAAA_,NVA*NOA*NOA*NVA,0.0);

  // Allocate and zero out Itermediates
  T * I_mnis_A   = this->memManager_->template malloc<T>(NB*NB*NB*NOA);
  T * I_mnij_AA  = this->memManager_->template malloc<T>(NB*NB*NOA*NOA);
  T * I_mbij_AAA = this->memManager_->template malloc<T>(NB*NOA*NOA*NVA);

  std::fill_n(I_mnis_A  ,NB*NB*NB*NOA,0.0);
  std::fill_n(I_mnij_AA ,NB*NB*NOA*NOA,0.0);
  std::fill_n(I_mbij_AAA,NB*NOA*NOA*NVA,0.0);

  // Conjugate for contraction if complex
  if(std::is_same<dcomplex,T>::value) {
    (*this->wfn_->moA()) = this->wfn_->moA()->conjugate();
  }

  // First Quarter (AAAA) transformation (mn|ls) -> (mn|is) (A)
  rank4w2Contract(3,this->wfn_->moA()->data(),NB,
    &this->wfn_->aointegrals()->aoERI_->storage()[0],NB,NB,NB,NB,
    I_mnis_A,NOA);

  // Conjugate back for contraction if complex
  if(std::is_same<dcomplex,T>::value) {
    (*this->wfn_->moA()) = this->wfn_->moA()->conjugate();
  }

  // Second Quarter (AAAA) transformation (mn|is) (A) -> (mn|ij) (AA)
  rank4w2Contract(4,this->wfn_->moA()->data(),NB,
    I_mnis_A,NB,NB,NOA,NB,
    I_mnij_AA,NOA);


  // Third Quarter (AAAA) transformation (mn|ij) (AA) -> (mb|ij) (AAA)
  rank4w2Contract(2,this->wfn_->moA()->data() + NB*NOA,NB,
    I_mnij_AA,NB,NB,NOA,NOA,
    I_mbij_AAA,NVA);

  // Conjugate for contraction if complex
  if(std::is_same<dcomplex,T>::value) {
    (*this->wfn_->moA()) = this->wfn_->moA()->conjugate();
  }

  // Fourth Quarter (AAAA) transformation (mb|ij) (AAA) -> (ab|ij) (AAAA)
  rank4w2Contract(1,this->wfn_->moA()->data() + NB*NOA,NB,
    I_mbij_AAA,NB,NVA,NOA,NOA,
    VVOOAAAA_,NVA);

  // Conjugate back for contraction if complex
  if(std::is_same<dcomplex,T>::value) {
    (*this->wfn_->moA()) = this->wfn_->moA()->conjugate();
  }

  // Free intermediate 1 and 3
  this->memManager_->free(I_mnis_A  ,NB*NB*NB*NOA);
  this->memManager_->free(I_mbij_AAA,NB*NOA*NOA*NVA);

  // If its closed shell, free intermediate 2
  if(!isOpenShell)
    this->memManager_->free(I_mnij_AA ,NB*NB*NOA*NOA);

  if(isOpenShell) {
    
    /// Evaluate VVOO (BBAA) MO Integrals

    // Allocate and Zero out VVOO (BBAA) storage
    VVOOBBAA_ = this->memManager_->template malloc<T>(NVB*NOA*NOA*NVB);
    std::fill_n(VVOOBBAA_,NVB*NOA*NOA*NVB,0.0);

    // Allocate and zero out Itermediates
    T * I_mbij_BAA = this->memManager_->template malloc<T>(NB*NOA*NOA*NVB);
    std::fill_n(I_mbij_BAA,NB*NOA*NOA*NVB,0.0);

    // Third Quarter (BBAA) transformation (mn|ij) (AA) -> (mb|ij) (BAA)
    rank4w2Contract(2,this->wfn_->moB()->data() + NB*NOB,NB,
      I_mnij_AA,NB,NB,NOA,NOA,
      I_mbij_BAA,NVB);

    // Conjugate for contraction if complex
    if(std::is_same<dcomplex,T>::value) {
      (*this->wfn_->moB()) = this->wfn_->moB()->conjugate();
    }
    
    // Fourth Quarter (BBAA) transformation (mb|ij) (BAA) -> (ab|ij) (BBAA)
    rank4w2Contract(1,this->wfn_->moB()->data() + NB*NOB,NB,
      I_mbij_BAA,NB,NVB,NOA,NOA,
      VVOOBBAA_,NVB);

    // Conjugate back for contraction if complex
    if(std::is_same<dcomplex,T>::value) {
      (*this->wfn_->moB()) = this->wfn_->moB()->conjugate();
    }


    // Free intermediates 2 and 3
    this->memManager_->free(I_mnij_AA ,NB*NB*NOA*NOA);
    this->memManager_->free(I_mbij_BAA,NB*NOA*NOA*NVB);

    /// Evaluate VVOO (BBBB) MO Integrals

    // Allocate and zero out VVOO (BBBB) storage
    VVOOBBBB_ = this->memManager_->template malloc<T>(NVB*NOB*NOB*NVB);
    std::fill_n(VVOOBBBB_,NVB*NOB*NOB*NVB,0.0);

    // Allocate and zero out Itermediates
    T * I_mnis_B   = this->memManager_->template malloc<T>(NB*NB*NB*NOB);
    T * I_mnij_BB  = this->memManager_->template malloc<T>(NB*NB*NOB*NOB);
    T * I_mbij_BBB = this->memManager_->template malloc<T>(NB*NOB*NOB*NVB);

    std::fill_n(I_mnis_B  ,NB*NB*NB*NOB,0.0);
    std::fill_n(I_mnij_BB ,NB*NB*NOB*NOB,0.0);
    std::fill_n(I_mbij_BBB,NB*NOB*NOB*NVB,0.0);

    // Conjugate for contraction if complex
    if(std::is_same<dcomplex,T>::value) {
      (*this->wfn_->moB()) = this->wfn_->moB()->conjugate();
    }

    // First Quarter (BBBB) transformation (mn|ls) -> (mn|is) (B)
    rank4w2Contract(3,this->wfn_->moB()->data(),NB,
      &this->wfn_->aointegrals()->aoERI_->storage()[0],NB,NB,NB,NB,
      I_mnis_B,NOB);

    // Conjugate back for contraction if complex
    if(std::is_same<dcomplex,T>::value) {
      (*this->wfn_->moB()) = this->wfn_->moB()->conjugate();
    }

    // Second Quarter (BBBB) transformation (mn|is) (B) -> (mn|ij) (BB)
    rank4w2Contract(4,this->wfn_->moB()->data(),NB,
      I_mnis_B,NB,NB,NOB,NB,
      I_mnij_BB,NOB);

    // Third Quarter (BBBB) transformation (mn|ij) (BB) -> (mb|ij) (BBB)
    rank4w2Contract(2,this->wfn_->moB()->data() + NB*NOB,NB,
      I_mnij_BB,NB,NB,NOB,NOB,
      I_mbij_BBB,NVB);

    // Conjugate for contraction if complex
    if(std::is_same<dcomplex,T>::value) {
      (*this->wfn_->moB()) = this->wfn_->moB()->conjugate();
    }

    // Fourth Quarter (BBBB) transformation (mb|ij) (BBB) -> (ab|ij) (BBBB)
    rank4w2Contract(1,this->wfn_->moB()->data() + NB*NOB,NB,
      I_mbij_BBB,NB,NVB,NOB,NOB,
      VVOOBBBB_,NVB);

    // Conjugate back for contraction if complex
    if(std::is_same<dcomplex,T>::value) {
      (*this->wfn_->moB()) = this->wfn_->moB()->conjugate();
    }

    this->memManager_->free(I_mnis_B  ,NB*NB*NB*NOB);
    this->memManager_->free(I_mbij_BBB,NB*NOB*NOB*NVB);

    VVOOAABB_ = this->memManager_->template malloc<T>(NVA*NOB*NOB*NVA);
    std::fill_n(VVOOAABB_,NVA*NOB*NOB*NVA,0.0);

    T * I_mbij_ABB = this->memManager_->template malloc<T>(NB*NOB*NOB*NVA);
    std::fill_n(I_mbij_ABB,NOB*NOB*NVA*NB,0.0);


    // Third Quarter (AABB) transformation (mn|ij) (BB) -> (mb|ij) (ABB)
    rank4w2Contract(2,this->wfn_->moA()->data() + NB*NOA,NB,
      I_mnij_BB,NB,NB,NOB,NOB,
      I_mbij_ABB,NVA);
    
    // Conjugate for contraction if complex
    if(std::is_same<dcomplex,T>::value) {
      (*this->wfn_->moA()) = this->wfn_->moA()->conjugate();
    }

    // Fourth Quarter (AABB) transformation (mb|ij) (ABB) -> (ab|ij) (AABB)
    rank4w2Contract(1,this->wfn_->moA()->data() + NB*NOA,NB,
      I_mbij_ABB,NB,NVA,NOB,NOB,
      VVOOAABB_,NVB);

    // Conjugate back for contraction if complex
    if(std::is_same<dcomplex,T>::value) {
      (*this->wfn_->moA()) = this->wfn_->moA()->conjugate();
    }

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

template <typename T> 
void MOIntegrals<T>::form2CVVOO() {

  int NB = this->wfn_->nBasis();
  int NO = this->wfn_->nO();
  int NV = this->wfn_->nV();

  VVOO_ = this->memManager_->template malloc<T>(NO*NV*NO*NV);
  std::fill_n(VVOO_,NO*NV*NO*NV,0.);

  T* tmp1 = this->memManager_->template malloc<T>(NO*NO*NB*NB);
  std::fill_n(tmp1,NO*NO*NB*NB,0.);

/*
  // N^8 Algorithm
  for(auto j = 0; j < NO; j++)
  for(auto i = 0; i < NO; i++)
  for(auto b = 0; b < NV; b++)
  for(auto a = 0; a < NV; a++)

  for(auto sg = 0; sg < 2*NB; sg+=2)
  for(auto lm = 0; lm < 2*NB; lm+=2)
  for(auto nu = 0; nu < 2*NB; nu+=2)
  for(auto mu = 0; mu < 2*NB; mu+=2){

    VVOO_[a + b*NV + i*NV*NV + j*NV*NO*NV] +=
      (
      std::conj((*wfn_->moA())(mu,a+NO)) * (*wfn_->moA())(nu,b+NO) *
      std::conj((*wfn_->moA())(lm,i)) * (*wfn_->moA())(sg,j) +

      std::conj((*wfn_->moA())(mu+1,a+NO)) * (*wfn_->moA())(nu+1,b+NO) *
      std::conj((*wfn_->moA())(lm,i))   * (*wfn_->moA())(sg,j) +

      std::conj((*wfn_->moA())(mu,a+NO))   * (*wfn_->moA())(nu,b+NO) *
      std::conj((*wfn_->moA())(lm+1,i)) * (*wfn_->moA())(sg+1,j) +

      std::conj((*wfn_->moA())(mu+1,a+NO)) * (*wfn_->moA())(nu+1,b+NO) *
      std::conj((*wfn_->moA())(lm+1,i)) * (*wfn_->moA())(sg+1,j) 
      ) * (*wfn_->aointegrals()->aoERI_)(mu/2,nu/2,lm/2,sg/2);
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

  this->wfn_->fileio()->out << "Begining First Half Transformation ABIJ"
    << endl;
  for(auto ij = 0; ij < NO*NO; ij+=nDo) {
    auto NDo = std::min(nDo, NO*NO-ij);
    for(auto iDo = 0, IJ = ij; iDo < NDo; iDo++, IJ++){ 
      int i = IJ % NO;
      int j = IJ / NO;
      scr2.noalias() = 
        wfn_->moA()->col(i).conjugate() * wfn_->moA()->col(j).transpose();
      scrMaps[iDo].noalias() = SAA + SBB;
    }

    TMap TMPMAP(tmp1+ij*NB*NB,NB*NB,NDo);
    TMPMAP.noalias() = aoERI * scr.block(0,0,NB*NB,NDo);
  }
 
  this->wfn_->fileio()->out << "Begining Second Half Transformation ABIJ"
    << endl;
  TMap TMPMAP(tmp1,NB*NB,NO*NO);
  TMap VVOOMAP(VVOO_,NV*NV,NO*NO);
  for(auto ab = 0; ab < NV*NV; ab+=nDo) {
    auto NDo = std::min(nDo, NV*NV - ab);
    for(auto iDo = 0, AB = ab; iDo < NDo; iDo++, AB++){ 
      int a = AB % NV;
      int b = AB / NV;
      scr2.noalias() = 
        wfn_->moA()->col(a+NO).conjugate() * wfn_->moA()->col(b+NO).transpose();
      scrMaps[iDo].noalias() = SAA + SBB;
    }

    VVOOMAP.block(ab,0,NDo,NO*NO).noalias() = 
      scr.block(0,0,NB*NB,NDo).transpose() * TMPMAP;
  }


  this->memManager_->free(tmp1,NO*NO*NB*NB);
  this->memManager_->free(scr.data(),nDo*NB*NB);
  this->memManager_->free(scr2.data(),4*NB*NB);
  this->haveMOVVOO_ = true;
}; // form2CVVOO

template <typename T>
void MOIntegrals<T>::formFullVVOO(){
  this->formVVOO();
  if(this->wfn_->nTCS() == 2) return;

  int NOA = this->wfn_->nOA();
  int NVA = this->wfn_->nVA();
  int NOB = this->wfn_->nOB();
  int NVB = this->wfn_->nVB();
  int NO = this->wfn_->nO();
  int NV = this->wfn_->nV();

  this->VVOO_ = this->memManager_->template malloc<T>(NV*NO*NV*NO);
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
} // formFullVVOO
