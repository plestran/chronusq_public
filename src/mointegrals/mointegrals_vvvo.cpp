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
#include <mointegrals.h>

namespace ChronusQ {
template<>
void MOIntegrals<double>::formVVVO(){
  if(this->haveMOVVVO_) return;

  bool is2C        = this->wfn_->nTCS() == 2;
  bool isOpenShell = !this->wfn_->isClosedShell;
  int NB           = this->wfn_->nBasis();
  int NOA = this->wfn_->nOA();
  int NVA = this->wfn_->nVA();
  int NOB = this->wfn_->nOB();
  int NVB = this->wfn_->nVB();


  /// Evaluate VVVO (AAAA) MO Integrals

  // Allocate and zero out VVVO (AAAA) storage
  VVVOAAAA_ = this->memManager_->malloc<double>(NVA*NVA*NVA*NOA);
  std::fill_n(VVVOAAAA_,NVA*NVA*NVA*NOA,0.0);

  // Allocate and zero out Itermediates
  double * I_1_A   = this->memManager_->malloc<double>(NB*NB*NB*NOA);
  double * I_2_AA  = this->memManager_->malloc<double>(NB*NB*NVA*NOA);
  double * I_3_AAA = this->memManager_->malloc<double>(NB*NVA*NVA*NOA);

  std::fill_n(I_1_A  ,NB*NB*NB*NOA,0.0);
  std::fill_n(I_2_AA ,NB*NB*NVA*NOA,0.0);
  std::fill_n(I_3_AAA,NB*NVA*NVA*NOA,0.0);

  // First Quarter (AAAA) transformation (mn|ls) -> (mn|li) (A)
  rank4w2Contract(4,this->wfn_->moA()->data(),NB,
    &this->wfn_->aointegrals()->aoERI_->storage()[0],NB,NB,NB,NB,
    I_1_A,NOA);

  // Second Quarter (AAAA) transformation (mn|li) (A) -> (mn|ci) (AA)
  rank4w2Contract(3,this->wfn_->moA()->data() + NOA*NB,NB,  
    I_1_A,NB,NB,NB,NOA,
    I_2_AA,NVA);

  // Third Quarter (AAAA) transformation (mn|ci) (AA) -> (mb|ci) (AAA)
  rank4w2Contract(2,this->wfn_->moA()->data() + NOA*NB,NB,  
    I_2_AA,NB,NB,NVA,NOA,
    I_3_AAA,NVA);

  // Fourth Quarter (AAAA) transformation (mb|ci) (AAA) -> (ab|ci) (AAAA)
  rank4w2Contract(1,this->wfn_->moA()->data() + NOA*NB,NB,  
    I_3_AAA,NB,NVA,NVA,NOA,
    VVVOAAAA_,NVA);

  this->memManager_->free(I_3_AAA,NB*NVA*NVA*NOA);
  this->memManager_->free(I_1_A  ,NB*NB*NB*NOA);

  if(isOpenShell) {
    // Allocate and zero out VVVO (BBAA) storage
    VVVOBBAA_ = this->memManager_->malloc<double>(NVB*NVB*NVA*NOA);
    std::fill_n(VVVOBBAA_,NVB*NVB*NVA*NOA,0.0);

    // Allocate and zero out Itermediates
    double * I_3_BAA = this->memManager_->malloc<double>(NB*NVB*NVA*NOA);

    std::fill_n(I_3_BAA,NB*NVB*NVA*NOA,0.0);

    // Third Quarter (AAAA) transformation (mn|ci) (AA) -> (mb|ci) (BAA)
    rank4w2Contract(2,this->wfn_->moB()->data() + NOB*NB,NB,  
      I_2_AA,NB,NB,NVA,NOA,
      I_3_BAA,NVB);
 
    // Fourth Quarter (AAAA) transformation (mb|ci) (BAA) -> (ab|ci) (BBAA)
    rank4w2Contract(1,this->wfn_->moB()->data() + NOB*NB,NB,  
      I_3_BAA,NB,NVB,NVA,NOA,
      VVVOBBAA_,NVB);

    this->memManager_->free(I_3_BAA,NB*NVB*NVA*NOA);

  }

  this->memManager_->free(I_2_AA ,NB*NB*NVA*NOA);

  if(isOpenShell) {
    /// Evaluate VVVO (BBBB) MO Integrals
 
    // Allocate and zero out VVVO (BBBB) storage
    VVVOBBBB_ = this->memManager_->malloc<double>(NVB*NVB*NVB*NOB);
    std::fill_n(VVVOBBBB_,NVB*NVB*NVB*NOB,0.0);
 
    // Allocate and zero out Itermediates
    double * I_1_B   = this->memManager_->malloc<double>(NB*NB*NB*NOB);
    double * I_2_BB  = this->memManager_->malloc<double>(NB*NB*NVB*NOB);
    double * I_3_BBB = this->memManager_->malloc<double>(NB*NVB*NVB*NOB);
 
    std::fill_n(I_1_B  ,NB*NB*NB*NOB,0.0);
    std::fill_n(I_2_BB ,NB*NB*NVB*NOB,0.0);
    std::fill_n(I_3_BBB,NB*NVB*NVB*NOB,0.0);
 
    // First Quarter (BBBB) transformation (mn|ls) -> (mn|li) (B)
    rank4w2Contract(4,this->wfn_->moB()->data(),NB,
      &this->wfn_->aointegrals()->aoERI_->storage()[0],NB,NB,NB,NB,
      I_1_B,NOB);
 
    // Second Quarter (BBBB) transformation (mn|li) (B) -> (mn|ci) (BB)
    rank4w2Contract(3,this->wfn_->moB()->data() + NOB*NB,NB,  
      I_1_B,NB,NB,NB,NOB,
      I_2_BB,NVB);
 
    // Third Quarter (BBBB) transformation (mn|ci) (BB) -> (mb|ci) (BBB)
    rank4w2Contract(2,this->wfn_->moB()->data() + NOB*NB,NB,  
      I_2_BB,NB,NB,NVB,NOB,
      I_3_BBB,NVB);
 
    // Fourth Quarter (BBBB) transformation (mb|ci) (BBB) -> (ab|ci) (BBBB)
    rank4w2Contract(1,this->wfn_->moB()->data() + NOB*NB,NB,  
      I_3_BBB,NB,NVB,NVB,NOB,
      VVVOBBBB_,NVB);
 
    this->memManager_->free(I_3_BBB,NB*NVB*NVB*NOB);
    this->memManager_->free(I_1_B  ,NB*NB*NB*NOB);


    // Allocate and zero out VVVO (AABB) storage
    VVVOAABB_ = this->memManager_->malloc<double>(NVA*NVA*NVB*NOB);
    std::fill_n(VVVOAABB_,NVA*NVA*NVB*NOB,0.0);

    // Allocate and zero out Itermediates
    double * I_3_ABB = this->memManager_->malloc<double>(NB*NVA*NVB*NOB);

    std::fill_n(I_3_ABB,NB*NVA*NVB*NOB,0.0);

    // Third Quarter (AABB) transformation (mn|ci) (BB) -> (mb|ci) (ABB)
    rank4w2Contract(2,this->wfn_->moA()->data() + NOA*NB,NB,  
      I_2_AA,NB,NB,NVB,NOB,
      I_3_ABB,NVA);
 
    // Fourth Quarter (AABB) transformation (mb|ci) (ABB) -> (ab|ci) (AABB)
    rank4w2Contract(1,this->wfn_->moA()->data() + NOA*NB,NB,  
      I_3_ABB,NB,NVA,NVB,NOB,
      VVVOAABB_,NVA);

    this->memManager_->free(I_3_ABB,NB*NVA*NVB*NOB);
    this->memManager_->free(I_2_BB ,NB*NB*NVB*NOB);
  }
}

template<>
void MOIntegrals<double>::formFullVVVO(){
  this->formVVVO();
  int NOA = this->wfn_->nOA();
  int NVA = this->wfn_->nVA();
  int NOB = this->wfn_->nOB();
  int NVB = this->wfn_->nVB();
  int NO = this->wfn_->nO();
  int NV = this->wfn_->nV();

  this->VVVO_ = this->memManager_->malloc<double>(NV*NV*NV*NO);
  std::fill_n(this->VVVO_,NO*NV*NV*NV,0.0);

  for(auto I = 0; I < NO; I+=2)
  for(auto C = 0; C < NV; C+=2) 
  for(auto B = 0; B < NV; B+=2) 
  for(auto A = 0; A < NV; A+=2){
    int a = A / 2;
    int b = B / 2;
    int c = C / 2;
    int i = I / 2;
    this->VVVO_[A + B*NV + C*NV*NV + I*NV*NV*NV] = 
      this->VVVOAAAA_[a + b*NVA + c*NVA*NVA + i*NVA*NVA*NVA];

    this->VVVO_[(A+1) + (B+1)*NV + C*NV*NV + I*NV*NV*NV] = 
      this->VVVOBBAA_[a + b*NVB + c*NVB*NVB + i*NVB*NVB*NVA];

    this->VVVO_[A + B*NV + (C+1)*NV*NV + (I+1)*NV*NV*NV] = 
      this->VVVOAABB_[a + b*NVA + c*NVA*NVA + i*NVA*NVA*NVB];

    this->VVVO_[(A+1) + (B+1)*NV + (C+1)*NV*NV + (I+1)*NV*NV*NV] = 
      this->VVVOBBBB_[a + b*NVB + c*NVB*NVB + i*NVB*NVB*NVB];
  }
} // formFullVVVO
};
