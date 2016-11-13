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
template<typename Op, typename T> 
void AOIntegrals::newTwoEContractDirect(
    const std::vector<std::reference_wrapper<Op>> &X,
    std::vector<std::reference_wrapper<Op>> &AX,
    std::vector<ERI_CONTRACTION_TYPE> &contractionList,
    std::vector<T>& scalingFactors){

  assert(X.size() == AX.size());
  assert(X.size() == contractionList.size());
  assert(X.size() == scalingFactors.size());

  size_t nMat = X.size();

#ifdef _OPENMP
  int nthreads = omp_get_max_threads();
#elif defined CQ_ENABLE_MPI
  int nthreads = getSize();
#else
  int nthreads = 1;
#endif

  // If the schwartz bound tensor hasn't been computed, compute it
  if(!this->haveSchwartz) this->computeSchwartz();
  // If the basis maps have not been computed, compute them
  if(!this->basisSet_->haveMapSh2Bf) this->basisSet_->makeMapSh2Bf(); 


  typedef Eigen::Matrix<typename Op::Scalar,Dynamic,Dynamic> TMatrix;
  typedef Eigen::Map<TMatrix> TMap;
  // Allocate space to hold thread specific contribution to
  std::vector<std::vector<TMap>> G;
  for(auto iMat = 0; iMat < nMat; iMat++){
    G.emplace_back();
  for(auto ithread = 0; ithread < nthreads; ithread++){
    G[iMat].emplace_back(
        this->memManager_->malloc<typename Op::Scalar>(this->nBasis_*this->nBasis_),
        this->nBasis_,this->nBasis_);
    G[iMat].back().setZero();
  };
  };
//  std::vector<std::vector<TMatrix>>(nMat,std::vector<TMatrix

  // each thread gets its own engine
  std::vector<libint2::Engine> engines(nthreads);

  // contruct engine for thread 0
  engines[0] = 
    libint2::Engine(libint2::Operator::coulomb,this->basisSet_->maxPrim(),
        this->basisSet_->maxL(),0);
  engines[0].set_precision(std::numeric_limits<double>::epsilon());
  
  // copy thread 0 engine to all others
  for(int i=1; i < nthreads; i++) engines[i] = engines[0];


  std::chrono::high_resolution_clock::time_point start,finish;

  // Allocate shBlkNorm storage
  this->basisSet_->shBlkNormAlpha = std::unique_ptr<RealMatrix>(
    new RealMatrix(this->basisSet_->nShell(),this->basisSet_->nShell())
  );

  std::vector<RealMatrix> shBlkNorms(nMat,
      RealMatrix::Zero(this->basisSet_->nShell(),this->basisSet_->nShell())
  );

  for(auto i = 0; i < nMat; i++)
    this->basisSet_->computeShBlkNorm(X[i].get(),shBlkNorms[i]);


  auto efficient_twoe_contract = [&] (int thread_id) {
    libint2::Engine &engine = engines[thread_id];
    for(auto s1 = 0ul, s1234 = 0ul; s1 < this->basisSet_->nShell(); s1++){
      size_t bf1_s = this->basisSet_->mapSh2Bf(s1);
      size_t n1    = this->basisSet_->shells(s1).size();
    for(auto s2 = 0ul; s2 <= s1; s2++){
      size_t bf2_s = this->basisSet_->mapSh2Bf(s2);
      size_t n2    = this->basisSet_->shells(s2).size();
    for(auto s3 = 0ul; s3 <= s1; s3++){
      size_t bf3_s = this->basisSet_->mapSh2Bf(s3);
      size_t n3    = this->basisSet_->shells(s3).size();
      size_t s4_max = (s1 == s3) ? s2 : s3;
    for(auto s4 = 0ul; s4 <= s4_max; s4++, s1234++){
      if(s1234 % nthreads != thread_id) continue;
      size_t bf4_s = this->basisSet_->mapSh2Bf(s4);
      size_t n4    = this->basisSet_->shells(s4).size();

      // Screening (Schwartz and Shell Block Norms)
      double shMax = 
        std::max(shBlkNorms[0](s1,s4),
        std::max(shBlkNorms[0](s2,s4),
        std::max(shBlkNorms[0](s3,s4),
        std::max(shBlkNorms[0](s1,s3),
        std::max(shBlkNorms[0](s2,s3),
                 shBlkNorms[0](s1,s2))))))
          * (*this->schwartz_)(s1,s2) * (*this->schwartz_)(s3,s4);
      for(auto iMat = 1ul; iMat < nMat; iMat++){
        shMax = std::max(shMax,
          std::max(shBlkNorms[iMat](s1,s4),
          std::max(shBlkNorms[iMat](s2,s4),
          std::max(shBlkNorms[iMat](s3,s4),
          std::max(shBlkNorms[iMat](s1,s3),
          std::max(shBlkNorms[iMat](s2,s3),
                   shBlkNorms[iMat](s1,s2))))))
            * (*this->schwartz_)(s1,s2) * (*this->schwartz_)(s3,s4));
      };

      if(shMax < this->thresholdSchwartz_ ) continue;

      double s12_deg = (s1 == s2) ? 1.0 : 2.0;
      double s34_deg = (s3 == s4) ? 1.0 : 2.0;
      double s12_34_deg = (s1 == s3) ? (s2 == s4 ? 1.0 : 2.0) : 2.0;
      double s1234_deg = s12_deg * s34_deg * s12_34_deg;

      const double * buff = engine.compute(
        this->basisSet_->shells(s1),
        this->basisSet_->shells(s2),
        this->basisSet_->shells(s3),
        this->basisSet_->shells(s4)
      );

      for(auto iMat = 0; iMat < nMat; iMat++) {
        for(auto i = 0ul, bf1 = bf1_s, ijkl = 0ul; i < n1; ++i, ++bf1         )
        for(auto j = 0ul, bf2 = bf2_s            ; j < n2; ++j, ++bf2         )
        for(auto k = 0ul, bf3 = bf3_s            ; k < n3; ++k, ++bf3         )
        for(auto l = 0ul, bf4 = bf4_s            ; l < n4; ++l, ++bf4, ++ijkl ){
          if(scalingFactors[iMat] == 0.0) continue;
          double v = buff[ijkl] * s1234_deg;

          if(contractionList[iMat] == COULOMB) {
            G[iMat][thread_id](bf1,bf2) += reinterpret_cast<double(&)[2]>(X[iMat].get()(bf4,bf3))[0] * v;
            G[iMat][thread_id](bf3,bf4) += reinterpret_cast<double(&)[2]>(X[iMat].get()(bf2,bf1))[0] * v;
            G[iMat][thread_id](bf2,bf1) += reinterpret_cast<double(&)[2]>(X[iMat].get()(bf3,bf4))[0] * v;
            G[iMat][thread_id](bf4,bf3) += reinterpret_cast<double(&)[2]>(X[iMat].get()(bf1,bf2))[0] * v;
          } else if(contractionList[iMat] == EXCHANGE) {
            G[iMat][thread_id](bf1,bf3) += 0.5 * X[iMat].get()(bf2,bf4) * v;
            G[iMat][thread_id](bf2,bf4) += 0.5 * X[iMat].get()(bf1,bf3) * v;
            G[iMat][thread_id](bf1,bf4) += 0.5 * X[iMat].get()(bf2,bf3) * v;
            G[iMat][thread_id](bf2,bf3) += 0.5 * X[iMat].get()(bf1,bf4) * v;
  
            G[iMat][thread_id](bf3,bf1) += 0.5 * X[iMat].get()(bf4,bf2) * v;
            G[iMat][thread_id](bf4,bf2) += 0.5 * X[iMat].get()(bf3,bf1) * v;
            G[iMat][thread_id](bf4,bf1) += 0.5 * X[iMat].get()(bf3,bf2) * v;
            G[iMat][thread_id](bf3,bf2) += 0.5 * X[iMat].get()(bf4,bf1) * v;
          }

        } // loop over BFs
      }; // loop over iMat

    }; // s4
    }; // s3
    }; // s2
    }; // s1
  }; // efficient_twoe_contract

  // Perform contraction (possibly in parallel)
  #pragma omp parallel
  {
    int thread_id = omp_get_thread_num();
    efficient_twoe_contract(thread_id);
  }

//  for(auto iMat = 0; iMat < nMat; iMat++) AX[iMat].get().setZero();

  // Reduce
  for(auto iMat = 0; iMat < nMat    ; iMat++)
  for(auto iTh  = 0; iTh  < nthreads; iTh++ ){
    AX[iMat].get() += 0.5 * scalingFactors[iMat] * G[iMat][iTh];
    this->memManager_->free(G[iMat][iTh].data(),G[iMat][iTh].size());
  }
};


template<typename Op, typename T> 
void AOIntegrals::newTwoEContractIncore(
    const std::vector<std::reference_wrapper<Op>> &X,
    std::vector<std::reference_wrapper<Op>> &AX,
    std::vector<ERI_CONTRACTION_TYPE> &contractionList,
    std::vector<T>& scalingFactors){

  typedef typename Op::Scalar Field;

  assert(X.size() == AX.size());
  assert(X.size() == contractionList.size());
  assert(X.size() == scalingFactors.size());

  size_t nMat = X.size();
  if(!this->haveAOTwoE) this->computeAOTwoE();

  RealTensor2d XTensor(X[0].get().cols(),X[0].get().rows());
  RealTensor2d AXTensor(X[0].get().cols(),X[0].get().rows());
  RealTensor2d IXTensor, IAXTensor;
  if(typeid(Field).hash_code() == typeid(dcomplex).hash_code()){
    IXTensor  = RealTensor2d(X[0].get().cols(),X[0].get().rows());
    IAXTensor = RealTensor2d(X[0].get().cols(),X[0].get().rows());
  }
  enum{i,j,k,l}; 

//for(auto iMat = 0; iMat < nMat; iMat++)
//  AX[iMat].get().setZero();

  for(auto iMat = 0; iMat < nMat; iMat++) {
    if(scalingFactors[iMat] == 0.0) continue;
    XTensor.fill(0.0);
    AXTensor.fill(0.0);

    // Copy Matrix over
    for(auto I = 0; I < X[iMat].get().size(); I++) {
      XTensor.storage()[I] = 
        reinterpret_cast<double(&)[2]>(X[iMat].get().data()[I])[0];
      if(typeid(Field).hash_code() == typeid(dcomplex).hash_code())
        IXTensor.storage()[I] = 
          reinterpret_cast<double(&)[2]>(X[iMat].get().data()[I])[1];
    }

    if(contractionList[iMat] == COULOMB) 
      contract(scalingFactors[iMat],*this->aoERI_,
        {i,j,k,l},XTensor,{l,k},0.0,AXTensor,{i,j});
    else if(contractionList[iMat] == EXCHANGE) 
      contract(scalingFactors[iMat],*this->aoERI_,
        {i,l,k,j},XTensor,{l,k},0.0,AXTensor,{i,j});

    if(typeid(Field).hash_code() == typeid(dcomplex).hash_code()) {
      if(contractionList[iMat] == EXCHANGE) 
        contract(scalingFactors[iMat],*this->aoERI_,
          {i,l,k,j},IXTensor,{l,k},0.0,IAXTensor,{i,j});
    }

    // Copy Tensor over
    for(auto I = 0; I < X[iMat].get().size(); I++){
      reinterpret_cast<double(&)[2]>(AX[iMat].get().data()[I])[0] += 
        2.0 * AXTensor.storage()[I];
      if(typeid(Field).hash_code() == typeid(dcomplex).hash_code())
        reinterpret_cast<double(&)[2]>(AX[iMat].get().data()[I])[1] += 
          2.0 * IAXTensor.storage()[I];
    }

  };
};
