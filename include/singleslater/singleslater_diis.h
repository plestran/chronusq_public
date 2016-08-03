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

#include <singleslater/singleslater_olddiis.h>

template<typename T>
void SingleSlater<T>::CDIIS4(int NDIIS){
  int N = NDIIS + 1;
  int NRHS = 1;
  int INFO = -1;
  int NB = this->nBasis_; 
  int NBSq = NB*NB;
  char NORM = 'O';
  double RCOND;


//this->fileio_->out << std::setw(2) << " " <<
//  std::setw(4) << " " <<
//  "*** Performing CDIIS Extrapolation ***" << endl;
  TMap B(this->memManager_->template malloc<T>(N*N),N,N);
  T   * coef   = this->memManager_->template malloc<T>(N);
  int * iPiv   = this->memManager_->template malloc<int>(N);
  int * iWORK_ = this->memManager_->template malloc<int>(N);


  for(auto j = 0; j < NDIIS; j++){
    // Scalar Part 
    this->readDIIS(this->EScalarDIIS_,j,this->NBSqScratch_->data());
    for(auto k = 0; k <= j; k++){
      this->readDIIS(this->EScalarDIIS_,k,this->NBSqScratch2_->data());

      B(j,k) = this->NBSqScratch_->frobInner(*this->NBSqScratch2_);
    }

    // Vector Part
    if(this->nTCS_ == 2 || !this->isClosedShell) {
      // Mz
      this->readDIIS(this->EMzDIIS_,j,this->NBSqScratch_->data());
      for(auto k = 0; k <= j; k++){
        this->readDIIS(this->EMzDIIS_,k,this->NBSqScratch2_->data());
     
        B(j,k) += this->NBSqScratch_->frobInner(*this->NBSqScratch2_);
      }
    } // has vector part
  } // Loop over errors

  B = B.template selfadjointView<Lower>();

  for (auto l=0;l < N-1;l++){
     B(N-1,l)=-1.0;
     B(l,N-1)=-1.0;
  }

  B(N-1,N-1)=0;
  for(auto k = 0; k < N;k++) coef[k] = 0.0; 
  coef[N-1]=-1.0;

  TMap COEFF(coef,N,1);
  //prettyPrint(this->fileio_->out,B,"CDIIS B Metric",20);
  //prettyPrint(this->fileio_->out,COEFF,"CDIIS RHS");
  
  double ANORM = B.template lpNorm<1>();

  int LWORK  = 5*N;
  T   *WORK  = this->memManager_->template malloc<T>(LWORK);

  TMatrix Bp(B);
  if(typeid(T).hash_code() == typeid(dcomplex).hash_code()){
    int LRWORK = 3*N;
    double   *RWORK = this->memManager_->
      template malloc<double>(LRWORK);

    // Linear Solve
    zgesv_(&N,&NRHS,reinterpret_cast<dcomplex*>(B.data()),&N,
        iPiv,reinterpret_cast<dcomplex*>(coef),&N,&INFO);

    // Obtain condition number from LU given from ZGESV
    zgecon_(&NORM,&N,reinterpret_cast<dcomplex*>(B.data()),&N,
        &ANORM,&RCOND,reinterpret_cast<dcomplex*>(WORK),RWORK,
        &INFO);

    this->memManager_->free(RWORK,LRWORK);
  } else {
    // Linear Solve
    dgesv_(&N,&NRHS,reinterpret_cast<double*>(B.data()),&N,
        iPiv,reinterpret_cast<double*>(coef),&N,&INFO);

    // Obtain condition number from LU given from DGESV
    dgecon_(&NORM,&N,reinterpret_cast<double*>(B.data()),&N,
        &ANORM,&RCOND,reinterpret_cast<double*>(WORK),iWORK_,
        &INFO);
  }
  this->memManager_->template free(WORK,LWORK);
  //prettyPrint(this->fileio_->out,COEFF,"CDIIS SOULTION");
  //prettyPrint(this->fileio_->out,Bp*COEFF,"SOLN");

  
/*
  int NUSE = N;
  if(std::abs(RCOND)<std::numeric_limits<double>::epsilon()) {
    this->fileio_->out << std::setw(2) << " " <<
      std::setw(4) << " " <<
      "*** XGESV Inversion Failed (RCOND = " 
      << std::abs(RCOND) << ") ***" << endl;

    std::fill(coef,coef+N,T(0.0));
    coef[N-1]=-1.0;
    T* S = this->memManager_->template malloc<T>(N);
    int RANK; 
    RCOND = -1;
    dgelss_(&N,&N,&NRHS,reinterpret_cast<double*>(Bp.data()),
      &N,reinterpret_cast<double*>(coef),&N,
      reinterpret_cast<double*>(S),&RCOND,&RANK,
      reinterpret_cast<double*>(WORK),&LWORK,&INFO);
      
    this->memManager_->free(S,N);
    cout << "NEW RANK" << RANK << endl;
    prettyPrint(this->fileio_->out,COEFF,"New CDIIS SOULTION");
  }
*/

    if(this->nTCS_ == 1 && this->isClosedShell){
      this->fockA_->setZero();
      this->PTA_->setZero();
    } else {
      this->fockScalar_->setZero();
      this->fockMz_->setZero();
      this->PTScalar_->setZero();
      this->PTMz_->setZero();
      if(this->nTCS_ == 2) {
        this->fockMy_->setZero();
        this->fockMx_->setZero();
        this->PTMy_->setZero();
        this->PTMx_->setZero();
      }
    }

    // Extrapolate full Fock matrix
    for(auto j = 0; j < NDIIS; j++) {
      this->readDIIS(this->FScalarDIIS_,j,this->NBSqScratch_->data());

      if(this->nTCS_ == 1 && this->isClosedShell)
        // Scalar Part
        (*this->fockA_) += coef[j] * (*this->NBSqScratch_);
      else {
        // Scalar Part
        this->fockScalar_->noalias() += coef[j] * (*this->NBSqScratch_);

        // Mz
        this->readDIIS(this->FMzDIIS_,j,this->NBSqScratch_->data());
        this->fockMz_->noalias() += coef[j] * (*this->NBSqScratch_);
      } // both vector and scalar
    } // Fock loop

    // Extrapolate G[P] for energy evaluation
    for(auto j = 0; j < NDIIS; j++) {
      this->readDIIS(this->PTScalarDIIS_,j,this->NBSqScratch_->data());

      if(this->nTCS_ == 1 && this->isClosedShell)
        // Scalar Part
        (*this->PTA_) += coef[j] * (*this->NBSqScratch_);
      else {
        // Scalar Part
        this->PTScalar_->noalias() += coef[j] * (*this->NBSqScratch_);

        // Mz
        this->readDIIS(this->PTMzDIIS_,j,this->NBSqScratch_->data());
        this->PTMz_->noalias() += coef[j] * (*this->NBSqScratch_);
      } // both vector and scalar
    } // PT loop
/*
  if(this->nTCS_ == 1 && this->isClosedShell){
    this->PTA_->real() = this->fockA_->real() - (*this->aointegrals_->coreH_);
    if(this->isDFT)
      for(auto j = 0; j < NDIIS; j++){
        this->readDIIS(this->VxcScalarDIIS_,j,this->NBSqScratch_->data());
        (*this->PTA_) -= coef[j] * (*this->NBSqScratch_);
      }
  } else {
    this->PTScalar_->real() = 
      this->fockScalar_->real() - (*this->aointegrals_->coreH_);
  }
*/

  this->orthoFock3();

  this->memManager_->free(B.data(),N*N);
  this->memManager_->free(coef,N);
  this->memManager_->free(iPiv,N);
  this->memManager_->free(iWORK_,N);
}

template<typename F> inline F DIISComplexScale();
template<>
inline double DIISComplexScale<double>(){return -1.; }
template<>
inline dcomplex DIISComplexScale<dcomplex>(){return dcomplex(0,1); }

template<typename T>
void SingleSlater<T>::genDIISCom(int iter){
  hsize_t offset[] = {iter % (this->nDIISExtrap_),0,0};
  hsize_t stride[] = {1,1,1};
  hsize_t block[]  = {1,1,1};
 
  hsize_t subDim[] = {1,this->nBasis_,this->nBasis_};
  hsize_t count[]  = {1,this->nBasis_,this->nBasis_};

  H5::DataSpace E = this->EScalarDIIS_->getSpace();
  E.selectHyperslab(H5S_SELECT_SET,count,offset,stride,block);

  H5::DataSpace memSpace(3,subDim,NULL);

  // Scalar Part
  // E(S) = [F(S),D(S)] + [F(K),D(K)]
  if(this->nTCS_ == 1 && this->isClosedShell){
//  this->NBSqScratch_->noalias() = 
//    (*this->fockOrthoA_) * (*this->onePDMOrthoA_);
//  this->NBSqScratch_->noalias() -=
//    (*this->onePDMOrthoA_) * (*this->fockOrthoA_);
    this->NBSqScratch_->noalias() = 
      (*this->fockA_) * (*this->onePDMA_) * (*this->aointegrals_->overlap_);
    this->NBSqScratch_->noalias() -=
      (*this->aointegrals_->overlap_) * (*this->onePDMA_) * (*this->fockA_);
  } else {
    this->NBSqScratch_->noalias() = 
      (*this->fockOrthoScalar_) * (*this->onePDMOrthoScalar_);
    this->NBSqScratch_->noalias() -=
      (*this->onePDMOrthoScalar_) * (*this->fockOrthoScalar_);

    this->NBSqScratch_->noalias() += 
      (*this->fockOrthoMz_) * (*this->onePDMOrthoMz_);
    this->NBSqScratch_->noalias() -=
      (*this->onePDMOrthoMz_) * (*this->fockOrthoMz_);

    if(this->nTCS_ == 2){
      this->NBSqScratch_->noalias() += 
        (*this->fockOrthoMy_) * (*this->onePDMOrthoMy_);
      this->NBSqScratch_->noalias() -=
        (*this->onePDMOrthoMy_) * (*this->fockOrthoMy_);

      this->NBSqScratch_->noalias() += 
        (*this->fockOrthoMx_) * (*this->onePDMOrthoMx_);
      this->NBSqScratch_->noalias() -=
        (*this->onePDMOrthoMx_) * (*this->fockOrthoMx_);
    }
  }

//prettyPrint(cout,*this->NBSqScratch_,"E1 FOR DIIS " + std::to_string(iter));
//this->aointegrals_->Ortho1Trans(*this->NBSqScratch_,*this->NBSqScratch2_);
//this->EScalarDIIS_->write(this->NBSqScratch2_->data(),H5PredType<T>(),
//  memSpace,E);
  this->EScalarDIIS_->write(this->NBSqScratch_->data(),H5PredType<T>(),
    memSpace,E);

  // Magnetization part
  // E(K) = [F(S),D(K)] + [F(K),D(S)] + i*([F(K+1),D(K+2)] - [F(K+2),D(K+1)])
  //
  // K = 1 (X)
  // K = 2 (Y)
  // K = 3 (Z)
  //
  // F(4) = F(X)
  // F(5) = F(Y)
  if(this->nTCS_ == 2 || !this->isClosedShell){
    // Mz Part
      
    // [F(S),D(Z)]
    this->NBSqScratch_->noalias() = 
      (*this->fockOrthoScalar_) * (*this->onePDMOrthoMz_);
    this->NBSqScratch_->noalias() -=
      (*this->onePDMOrthoMz_) * (*this->fockOrthoScalar_);

    // [F(Z),D(S)]
    this->NBSqScratch_->noalias() += 
      (*this->fockOrthoMz_) * (*this->onePDMOrthoScalar_);
    this->NBSqScratch_->noalias() -=
      (*this->onePDMOrthoScalar_) * (*this->fockOrthoMz_);
 
    if(this->nTCS_ == 2){
      // [F(X),D(Y)]
      this->NBSqScratch_->noalias() += DIISComplexScale<T>() *
        (*this->fockOrthoMx_) * (*this->onePDMOrthoMy_);
      this->NBSqScratch_->noalias() -= DIISComplexScale<T>() *
        (*this->onePDMOrthoMy_) * (*this->fockOrthoMx_);
     
      // [F(Y),D(X)]
      this->NBSqScratch_->noalias() -= DIISComplexScale<T>() *
        (*this->fockOrthoMy_) * (*this->onePDMOrthoMx_);
      this->NBSqScratch_->noalias() += DIISComplexScale<T>() *
        (*this->onePDMOrthoMx_) * (*this->fockOrthoMy_);
    }
      

    this->EMzDIIS_->write(this->NBSqScratch_->data(),H5PredType<T>(),
      memSpace,E);

    if(this->nTCS_ == 2){
      // Mx Part
        
      // [F(S),D(X)]
      this->NBSqScratch_->noalias() = 
        (*this->fockOrthoScalar_) * (*this->onePDMOrthoMx_);
      this->NBSqScratch_->noalias() -=
        (*this->onePDMOrthoMx_) * (*this->fockOrthoScalar_);

      // [F(X),D(S)]
      this->NBSqScratch_->noalias() += 
        (*this->fockOrthoMx_) * (*this->onePDMOrthoScalar_);
      this->NBSqScratch_->noalias() -=
        (*this->onePDMOrthoScalar_) * (*this->fockOrthoMx_);

      // [F(Y),D(Z)]
      this->NBSqScratch_->noalias() += DIISComplexScale<T>() *
        (*this->fockOrthoMy_) * (*this->onePDMOrthoMz_);
      this->NBSqScratch_->noalias() -= DIISComplexScale<T>() *
        (*this->onePDMOrthoMz_) * (*this->fockOrthoMy_);
     
      // [F(Z),D(Y)]
      this->NBSqScratch_->noalias() -= DIISComplexScale<T>() *
        (*this->fockOrthoMz_) * (*this->onePDMOrthoMy_);
      this->NBSqScratch_->noalias() += DIISComplexScale<T>() *
        (*this->onePDMOrthoMy_) * (*this->fockOrthoMz_);

      this->EMxDIIS_->write(this->NBSqScratch_->data(),H5PredType<T>(),
        memSpace,E);

      // My Part
        
      // [F(S),D(Y)]
      this->NBSqScratch_->noalias() = 
        (*this->fockOrthoScalar_) * (*this->onePDMOrthoMy_);
      this->NBSqScratch_->noalias() -=
        (*this->onePDMOrthoMy_) * (*this->fockOrthoScalar_);

      // [F(Y),D(S)]
      this->NBSqScratch_->noalias() += 
        (*this->fockOrthoMy_) * (*this->onePDMOrthoScalar_);
      this->NBSqScratch_->noalias() -=
        (*this->onePDMOrthoScalar_) * (*this->fockOrthoMy_);

      // [F(Z),D(X)]
      this->NBSqScratch_->noalias() += DIISComplexScale<T>() *
        (*this->fockOrthoMz_) * (*this->onePDMOrthoMx_);
      this->NBSqScratch_->noalias() -= DIISComplexScale<T>() *
        (*this->onePDMOrthoMx_) * (*this->fockOrthoMz_);
     
      // [F(X),D(Z)]
      this->NBSqScratch_->noalias() -= DIISComplexScale<T>() *
        (*this->fockOrthoMx_) * (*this->onePDMOrthoMz_);
      this->NBSqScratch_->noalias() += DIISComplexScale<T>() *
        (*this->onePDMOrthoMz_) * (*this->fockOrthoMx_);

      this->EMyDIIS_->write(this->NBSqScratch_->data(),H5PredType<T>(),
        memSpace,E);
    }
  }

};

template<typename T>
void SingleSlater<T>::initDIISFiles(){
  std::vector<hsize_t> dims;
  dims.push_back(this->nDIISExtrap_);
  dims.push_back(this->nBasis_);
  dims.push_back(this->nBasis_);

  this->FScalarDIIS_ = 
    this->fileio_->createScratchPartition(H5PredType<T>(),
      "Fock (Scalar) For DIIS Extrapoloation",dims);
  this->DScalarDIIS_ = 
    this->fileio_->createScratchPartition(H5PredType<T>(),
      "Density (Scalar) For DIIS Extrapoloation",dims);
  this->EScalarDIIS_ = 
    this->fileio_->createScratchPartition(H5PredType<T>(),
      "Error Metric [F,D] (Scalar) For DIIS Extrapoloation",dims);
  this->PTScalarDIIS_ = 
    this->fileio_->createScratchPartition(H5PredType<T>(),
      "PT (Scalar) For DIIS Extrapoloation",dims);

  if(this->nTCS_ == 2 || !this->isClosedShell) {
    this->FMzDIIS_ = 
      this->fileio_->createScratchPartition(H5PredType<T>(),
        "Fock (Mz) For DIIS Extrapoloation",dims);
    this->DMzDIIS_ = 
      this->fileio_->createScratchPartition(H5PredType<T>(),
        "Density (Mz) For DIIS Extrapoloation",dims);
    this->EMzDIIS_ = 
      this->fileio_->createScratchPartition(H5PredType<T>(),
        "Error Metric [F,D] (Mz) For DIIS Extrapoloation",dims);
    this->PTMzDIIS_ = 
      this->fileio_->createScratchPartition(H5PredType<T>(),
        "PT (Mz) For DIIS Extrapoloation",dims);
  }

  if(this->nTCS_ == 2) {
    this->FMyDIIS_ = 
      this->fileio_->createScratchPartition(H5PredType<T>(),
        "Fock (My) For DIIS Extrapoloation",dims);
    this->DMyDIIS_ = 
      this->fileio_->createScratchPartition(H5PredType<T>(),
        "Density (My) For DIIS Extrapoloation",dims);
    this->EMyDIIS_ = 
      this->fileio_->createScratchPartition(H5PredType<T>(),
        "Error Metric [F,D] (My) For DIIS Extrapoloation",dims);
    this->PTMyDIIS_ = 
      this->fileio_->createScratchPartition(H5PredType<T>(),
        "PT (My) For DIIS Extrapoloation",dims);

    this->FMxDIIS_ = 
      this->fileio_->createScratchPartition(H5PredType<T>(),
        "Fock (Mx) For DIIS Extrapoloation",dims);
    this->DMxDIIS_ = 
      this->fileio_->createScratchPartition(H5PredType<T>(),
        "Density (Mx) For DIIS Extrapoloation",dims);
    this->EMxDIIS_ = 
      this->fileio_->createScratchPartition(H5PredType<T>(),
        "Error Metric [F,D] (Mx) For DIIS Extrapoloation",dims);
    this->PTMxDIIS_ = 
      this->fileio_->createScratchPartition(H5PredType<T>(),
        "PT (Mx) For DIIS Extrapoloation",dims);
  }

};

template<typename T>
void SingleSlater<T>::cpyDenDIIS(int iter){
  T* ScalarPtr;
  if(this->isClosedShell && this->nTCS_ == 1) 
    ScalarPtr = this->onePDMA_->data();
  else
    ScalarPtr = this->onePDMScalar_->data();

  this->writeDIIS(this->DScalarDIIS_,iter % this->nDIISExtrap_,ScalarPtr);
  if(this->nTCS_ == 2 || !this->isClosedShell){
    this->writeDIIS(this->DMzDIIS_,iter % this->nDIISExtrap_,
      this->onePDMMz_->data());
    if(this->nTCS_ == 2){
      this->writeDIIS(this->DMyDIIS_,iter % this->nDIISExtrap_,
        this->onePDMMy_->data());
      this->writeDIIS(this->DMxDIIS_,iter % this->nDIISExtrap_,
        this->onePDMMx_->data());
    }
  }
}

template<typename T>
void SingleSlater<T>::cpyFockDIIS(int iter){
  T* ScalarPtr;
  if(this->isClosedShell && this->nTCS_ == 1) 
    ScalarPtr = this->fockA_->data();
  else
    ScalarPtr = this->fockScalar_->data();

  this->writeDIIS(this->FScalarDIIS_,iter % this->nDIISExtrap_,ScalarPtr);

  if(this->isClosedShell && this->nTCS_ == 1) 
    ScalarPtr = this->PTA_->data();
  else
    ScalarPtr = this->PTScalar_->data();

  this->writeDIIS(this->PTScalarDIIS_,iter % this->nDIISExtrap_,
    ScalarPtr);

  if(this->nTCS_ == 2 || !this->isClosedShell){
    this->writeDIIS(this->FMzDIIS_,iter % this->nDIISExtrap_,
      this->fockMz_->data());
    if(this->nTCS_ == 2){
      this->writeDIIS(this->FMyDIIS_,iter % this->nDIISExtrap_,
        this->fockMy_->data());
      this->writeDIIS(this->FMxDIIS_,iter % this->nDIISExtrap_,
        this->fockMx_->data());
    }
  }
}

template <typename T>
void SingleSlater<T>::readDIIS(H5::DataSet *DIISF, int IDIISIter, T *data){
  hsize_t offset[] = {IDIISIter,0,0};
  

  hsize_t stride[] = {1,1,1};
  hsize_t block[]  = {1,1,1};
 
  hsize_t subDim[] = {1,this->nBasis_,this->nBasis_};
  hsize_t count[]  = {1,this->nBasis_,this->nBasis_};

  H5::DataSpace DATA = DIISF->getSpace();
  H5::DataSpace memspace(3,subDim,NULL);

  DATA.selectHyperslab(H5S_SELECT_SET,count,offset,stride,block);

  DIISF->read(data,H5PredType<T>(),memspace,DATA);

}

template <typename T>
void SingleSlater<T>::writeDIIS(H5::DataSet *DIISF, int IDIISIter, T *data){
  hsize_t offset[] = {IDIISIter,0,0};
  

  hsize_t stride[] = {1,1,1};
  hsize_t block[]  = {1,1,1};
 
  hsize_t subDim[] = {1,this->nBasis_,this->nBasis_};
  hsize_t count[]  = {1,this->nBasis_,this->nBasis_};

  H5::DataSpace DATA = DIISF->getSpace();
  H5::DataSpace memspace(3,subDim,NULL);

  DATA.selectHyperslab(H5S_SELECT_SET,count,offset,stride,block);

  DIISF->write(data,H5PredType<T>(),memspace,DATA);

}
