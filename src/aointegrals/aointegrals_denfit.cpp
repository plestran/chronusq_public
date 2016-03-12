/* 
 *  The Chronus Quantum (ChronusQ) software package is high-performace 
 *  computational chemistry software with a strong emphasis on explicitly 
 *  time-dependent and post-SCF quantum mechanical methods.
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
#include <aointegrals.h>
using ChronusQ::AOIntegrals;
void AOIntegrals::computeAORII(){
  if(!this->haveSchwartz) this->computeSchwartz();
  if(!this->haveRIS)      this->computeAORIS(); 

#ifdef _OPENMP
  int nthreads = omp_get_max_threads();
#else
  int nthreads = 1;
#endif

  std::vector<coulombEngine> engines(nthreads);
  engines[0] = coulombEngine(
    std::max(this->basisSet_->maxPrim(),this->DFbasisSet_->maxPrim()),
    std::max(this->basisSet_->maxL(),this->DFbasisSet_->maxL()),0);
  engines[0].set_precision(std::numeric_limits<double>::epsilon());

  for(int i=1; i<nthreads; i++) engines[i] = engines[0];
  if(!this->basisSet_->haveMapSh2Bf) this->basisSet_->makeMapSh2Bf(1); 
  if(!this->DFbasisSet_->haveMapSh2Bf) this->DFbasisSet_->makeMapSh2Bf(1); 

#ifdef _OPENMP
  #pragma omp parallel
#endif
  {
#ifdef _OPENMP
    int thread_id = omp_get_thread_num();
#else
    int thread_id = 0;
#endif
  for(int s1 = 0, s123=0; s1 < this->basisSet_->nShell(); s1++) {
    int bf1_s = this->basisSet_->mapSh2Bf(s1);
    int n1    = this->basisSet_->shells(s1).size();
    for(int s2 = 0; s2 < this->basisSet_->nShell(); s2++) {
      int bf2_s = this->basisSet_->mapSh2Bf(s2);
      int n2    = this->basisSet_->shells(s2).size();
      for(int dfs = 0; dfs < this->DFbasisSet_->nShell(); dfs++,s123++) {
        if(s123 % nthreads != thread_id) continue;
        int dfbf3_s = this->DFbasisSet_->mapSh2Bf(dfs);
        int dfn3    = this->DFbasisSet_->shells(dfs).size();

        // Schwartz and Density screening
        if((*this->schwartz_)(s1,s2) * (*this->aoRIS_)(dfs,dfs)
            < this->controls_->thresholdSchawrtz ) continue;
 
        const double* buff = engines[thread_id].compute(
          this->basisSet_->shells(s1),
          this->basisSet_->shells(s2),
          this->DFbasisSet_->shells(dfs),
          libint2::Shell::unit());

        auto lower = {bf1_s,bf2_s,dfbf3_s};
        auto upper = {bf1_s+n1,bf2_s+n2,dfbf3_s+dfn3};
        auto view  = btas::make_view(
          this->aoRII_->range().slice(lower,upper),
          this->aoRII_->storage());
        std::copy(buff,buff+n1*n2*dfn3,view.begin());
      }
    }
  }
  } // omp parallel scope
  this->haveRII = true;
}

void AOIntegrals::computeAORIS(){
#ifdef _OPENMP
  int nthreads = omp_get_max_threads();
#else
  int nthreads = 1;
#endif
  this->haveRIS= true;
  std::vector<coulombEngine> engines(nthreads);
  engines[0] = coulombEngine( this->DFbasisSet_->maxPrim(), 
                              this->DFbasisSet_->maxL(),0);
  engines[0].set_precision(std::numeric_limits<double>::epsilon());

  for(int i=1; i<nthreads; i++) engines[i] = engines[0];
  if(!this->DFbasisSet_->haveMapSh2Bf) this->DFbasisSet_->makeMapSh2Bf(1); 

  RealMap aoRISMap(&this->aoRIS_->storage()[0],
    this->DFbasisSet_->nBasis(),this->DFbasisSet_->nBasis());

#ifdef _OPENMP
  #pragma omp parallel
#endif
  {
#ifdef _OPENMP
    int thread_id = omp_get_thread_num();
#else
    int thread_id = 0;
#endif
  for(int s1 = 0, s12=0; s1 < this->DFbasisSet_->nShell(); s1++) {
    int bf1_s = this->DFbasisSet_->mapSh2Bf(s1);
    int n1    = this->DFbasisSet_->shells(s1).size();
    for(int s2 = 0; s2 < this->DFbasisSet_->nShell(); s2++,s12++) {
      int bf2_s = this->DFbasisSet_->mapSh2Bf(s2);
      int n2    = this->DFbasisSet_->shells(s2).size();
 
      if(s12 % nthreads != thread_id) continue;

      const double* buff = engines[thread_id].compute(
        this->DFbasisSet_->shells(s1),
        libint2::Shell::unit(),
        this->DFbasisSet_->shells(s2),
        libint2::Shell::unit());

      ConstRealMap buffMat(buff,n1,n2);
      aoRISMap.block(bf1_s,bf2_s,n1,n2) = buffMat; 
    }
  }
  } // omp parallel scope

  aoRISMap = aoRISMap.selfadjointView<Lower>(); // Symmetrize
  this->haveRIS = true;
}

void AOIntegrals::transformAORII(){
  if(!this->haveRIS) this->computeAORIS();
  if(!this->haveRII) this->computeAORII();

  RealMap S(&this->aoRIS_->storage()[0],this->DFbasisSet_->nBasis(),this->DFbasisSet_->nBasis());
  RealMatrix ShalfMat = S.pow(-0.5);
  RealTensor2d Shalf(this->DFbasisSet_->nBasis(),this->DFbasisSet_->nBasis());
  for(auto i = 0; i < Shalf.size(); i++)
    Shalf.storage()[i] = ShalfMat.data()[i];
  RealTensor3d A0(this->nBasis_,this->nBasis_,this->DFbasisSet_->nBasis());
  enum{i,j,k,l,X,Y};
  contract(1.0,*this->aoRII_,{i,j,X},Shalf,{X,Y},0.0,A0,{i,j,Y});
  *this->aoRII_ = A0;
  this->haveTRII = true;
}

void AOIntegrals::compareRI(){
  if(!this->haveRIS)  this->computeAORIS();
  if(!this->haveRII)  this->computeAORII();
  if(!this->haveTRII) this->transformAORII();


  RealTensor4d fakeERI(this->nBasis_,this->nBasis_,this->nBasis_,this->nBasis_);
  RealTensor4d diff(this->nBasis_,this->nBasis_,this->nBasis_,this->nBasis_);

  enum{i,j,k,l,X,Y};
  contract(1.0,*this->aoRII_,{i,j,X},*this->aoRII_,{k,l,X},0.0,fakeERI,{i,j,k,l});

  diff = fakeERI - (*this->aoERI_);
  auto sum = 0.0;
  for(auto x : diff) sum += std::abs(x);
  cout << std::scientific << sum / diff.size() << endl;
}


