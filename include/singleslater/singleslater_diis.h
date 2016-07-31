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

template<typename T>
void SingleSlater<T>::CDIIS2(){
  int N = this->lenCoeff_;
  int NRHS = 1;
  int INFO = -1;
  int NB = this->nBasis_ * this->nTCS_; 
  int NBSq = NB*NB;
  char NORM = 'O';
  double RCOND;

  this->fileio_->out << std::setw(2) << " " <<
    std::setw(4) << " " <<
    "*** Performing CDIIS Extrapolation ***" << endl;
  TMap B(this->memManager_->template malloc<T>(N*N),N,N);
  T   * coef   = this->memManager_->template malloc<T>(N);
  int * iPiv   = this->memManager_->template malloc<int>(N);
  int * iWORK_ = this->memManager_->template malloc<int>(N);

/*
  for(auto j = 0; j < (N-1); j++)
  for(auto k = 0; k <= j          ; k++){
    TMap EJA(this->ErrorAlphaMem_ + (j%(N-1))*NBSq,NB,NB);
    TMap EKA(this->ErrorAlphaMem_ + (k%(N-1))*NBSq,NB,NB);
//    B(j,k) = -EJA.frobInner(EKA);
    B(j,k) = EJA.frobInner(EKA);
    if(!this->isClosedShell && this->nTCS_ == 1){
      TMap EJB(this->ErrorBetaMem_ + (j%(N-1))*NBSq,NB,NB);
      TMap EKB(this->ErrorBetaMem_ + (k%(N-1))*NBSq,NB,NB);
//      B(j,k) += -EJB.frobInner(EKB);
      B(j,k) += EJB.frobInner(EKB);
    }
    B(k,j) = B(j,k);
  }
*/

  hsize_t stride[] = {1,1,1};
  hsize_t block[]  = {1,1,1};
 
  hsize_t subDim[] = {1,this->nBasis_,this->nBasis_};
  hsize_t count[]  = {1,this->nBasis_,this->nBasis_};

  H5::DataSpace EJ = this->EScalarDIIS_->getSpace();
  H5::DataSpace EK = this->EScalarDIIS_->getSpace();
  H5::DataSpace memSpace(3,subDim,NULL);

  for(auto j = 0; j < (N-1); j++){
    hsize_t offset_j[] = {j,0,0};
    EJ.selectHyperslab(H5S_SELECT_SET,count,offset_j,stride,block);

    // Scalar Part
    this->EScalarDIIS_->read(this->NBSqScratch_->data(),H5PredType<T>(),
      memSpace,EJ);

    for(auto k = 0; k <= j; k++){
      hsize_t offset_k[] = {k,0,0};
      EK.selectHyperslab(H5S_SELECT_SET,count,offset_j,stride,block);
      this->EScalarDIIS_->read(this->NBSqScratch2_->data(),H5PredType<T>(),
        memSpace,EK);

      B(j,k) = -this->NBSqScratch_->frobInner(*this->NBSqScratch2_);
    }

    // Vector Part

    if(this->nTCS_ == 2 || !this->isClosedShell) {
      // Mz
      this->EMzDIIS_->read(this->NBSqScratch_->data(),H5PredType<T>(),
        memSpace,EJ);
      for(auto k = 0; k <= j; k++){
        hsize_t offset_k[] = {k,0,0};
        EK.selectHyperslab(H5S_SELECT_SET,count,offset_j,stride,block);
        this->EMzDIIS_->read(this->NBSqScratch2_->data(),H5PredType<T>(),
          memSpace,EK);
     
        B(j,k) += -this->NBSqScratch_->frobInner(*this->NBSqScratch2_);
      }

      if(this->nTCS_ == 2){
        // Mx
        this->EMxDIIS_->read(this->NBSqScratch_->data(),H5PredType<T>(),
          memSpace,EJ);
        for(auto k = 0; k <= j; k++){
          hsize_t offset_k[] = {k,0,0};
          EK.selectHyperslab(H5S_SELECT_SET,count,offset_j,stride,block);
          this->EMxDIIS_->read(this->NBSqScratch2_->data(),H5PredType<T>(),
            memSpace,EK);
     
          B(j,k) += -this->NBSqScratch_->frobInner(*this->NBSqScratch2_);
        }

        // My
        this->EMyDIIS_->read(this->NBSqScratch_->data(),H5PredType<T>(),
          memSpace,EJ);
        for(auto k = 0; k <= j; k++){
          hsize_t offset_k[] = {k,0,0};
          EK.selectHyperslab(H5S_SELECT_SET,count,offset_j,stride,block);
          this->EMyDIIS_->read(this->NBSqScratch2_->data(),H5PredType<T>(),
            memSpace,EK);
     
          B(j,k) += -this->NBSqScratch_->frobInner(*this->NBSqScratch2_);
        }

      } // ntcs = 2
    } // has vector part

  }

  B = B.template selfadjointView<Lower>();
  if(this->nTCS_ == 2 || !this->isClosedShell) B *= 0.25;

  for (auto l=0;l<N-1;l++){
     B(N-1,l)=-1.0;
     B(l,N-1)=-1.0;
  }

  B(N-1,N-1)=0;
  for(auto k = 0; k < N;k++) coef[k] = 0.0; 
  coef[N-1]=-1.0;

  TMap COEFF(coef,N,1);
  prettyPrint(this->fileio_->out,B,"CDIIS B Metric",20);
  prettyPrint(this->fileio_->out,COEFF,"CDIIS RHS");
  
  double ANORM = B.template lpNorm<1>();

  int LWORK  = 5*this->nDIISExtrap_;
  T   *WORK  = this->memManager_->template malloc<T>(LWORK);

  if(typeid(T).hash_code() == typeid(dcomplex).hash_code()){
    int LRWORK = 3*this->nDIISExtrap_;
    double   *RWORK = this->memManager_->template malloc<double>(LRWORK);
    zgesv_(&N,&NRHS,reinterpret_cast<dcomplex*>(B.data()),&N,iPiv,
        reinterpret_cast<dcomplex*>(coef),&N,&INFO);
    zgecon_(&NORM,&N,reinterpret_cast<dcomplex*>(B.data()),&N,&ANORM,&RCOND,
        reinterpret_cast<dcomplex*>(WORK),RWORK,&INFO);
    this->memManager_->free(RWORK,LRWORK);
  } else {
    dgesv_(&N,&NRHS,reinterpret_cast<double*>(B.data()),&N,iPiv,
        reinterpret_cast<double*>(coef),&N,&INFO);
    dgecon_(&NORM,&N,reinterpret_cast<double*>(B.data()),&N,&ANORM,&RCOND,
        reinterpret_cast<double*>(WORK),iWORK_,&INFO);
  }
  this->memManager_->template free(WORK,LWORK);


  if(std::abs(RCOND) > std::numeric_limits<double>::epsilon()) {
    this->fockA_->setZero();
    if(!this->isClosedShell && this->nTCS_ == 1) 
      this->fockB_->setZero();

/*
    for(auto j = 0; j < N-1; j++) {
      TMap FA(this->FADIIS_ + (j%(N-1))*NBSq,NB,NB);
      (*this->fockA_) += coef[j]*FA;
      if(!this->isClosedShell && this->nTCS_ == 1) {
        TMap FB(this->FBDIIS_ + (j%(N-1))*NBSq,NB,NB);
        (*this->fockB_) += coef[j]*FB;
      }
    }
*/
    for(auto j = 0; j < N-1; j++) {
      hsize_t offset_j[] = {j,0,0};
      EJ.selectHyperslab(H5S_SELECT_SET,count,offset_j,stride,block);
      this->FScalarDIIS_->read(this->NBSqScratch_->data(),H5PredType<T>(),
        memSpace,EJ);

      if(this->nTCS_ == 1 && this->isClosedShell)
        // Scalar Part
//      this->fockA_->noalias() += coef[j] * (*this->NBSqScratch_);
        this->fockOrthoA_->noalias() += coef[j] * (*this->NBSqScratch_);
      else {
        // Scalar Part
        this->fockScalar_->noalias() += coef[j] * (*this->NBSqScratch_);

        // Mz
        this->FMzDIIS_->read(this->NBSqScratch_->data(),H5PredType<T>(),
          memSpace,EJ);
        this->fockMz_->noalias() += coef[j] * (*this->NBSqScratch_);

        if(this->nTCS_ == 2){
          // My
          this->FMyDIIS_->read(this->NBSqScratch_->data(),H5PredType<T>(),
            memSpace,EJ);
          this->fockMy_->noalias() += coef[j] * (*this->NBSqScratch_);
          // Mx
          this->FMxDIIS_->read(this->NBSqScratch_->data(),H5PredType<T>(),
            memSpace,EJ);
          this->fockMx_->noalias() += coef[j] * (*this->NBSqScratch_);
        } // ntcs = 2
      } // both vector and scalar
    } // loop
  } else {
    this->fileio_->out << std::setw(2) << " " <<
      std::setw(4) << " " <<
      "*** CDIIS Extrapolation Failed (RCOND = " << std::abs(RCOND) << 
      ") ***" << endl;
  }

  // Gather matricies
  if(this->nTCS_ == 2 || !this->isClosedShell){
    std::vector<std::reference_wrapper<TMap>> toGather;
    toGather.emplace_back(*this->fockScalar_);
    toGather.emplace_back(*this->fockMz_);
    if(this->nTCS_ == 2){
      toGather.emplace_back(*this->fockMy_);
      toGather.emplace_back(*this->fockMx_);
      Quantum<T>::spinGather(*this->fockA_,toGather);
    } else 
      Quantum<T>::spinGather(*this->fockA_,*this->fockB_,toGather);
  }

  this->memManager_->free(B.data(),N*N);
  this->memManager_->free(coef,N);
  this->memManager_->free(iPiv,N);
  this->memManager_->free(iWORK_,N);

} // CDIIS

template<typename T>
void SingleSlater<T>::CDIIS(){
  int N = this->lenCoeff_;
  int NRHS = 1;
  int INFO = -1;
  int NB = this->nBasis_ * this->nTCS_; 
  int NBSq = NB*NB;
  char NORM = 'O';
  double RCOND;

  this->fileio_->out << std::setw(2) << " " <<
    std::setw(4) << " " <<
    "*** Performing CDIIS Extrapolation ***" << endl;
  TMap B(this->memManager_->template malloc<T>(N*N),N,N);
  T   * coef   = this->memManager_->template malloc<T>(N);
  int * iPiv   = this->memManager_->template malloc<int>(N);
  int * iWORK_ = this->memManager_->template malloc<int>(N);

  for(auto j = 0; j < (N-1); j++)
  for(auto k = 0; k <= j          ; k++){
    TMap EJA(this->ErrorAlphaMem_ + (j%(N-1))*NBSq,NB,NB);
    TMap EKA(this->ErrorAlphaMem_ + (k%(N-1))*NBSq,NB,NB);
//    B(j,k) = -EJA.frobInner(EKA);
    B(j,k) = EJA.frobInner(EKA);
    if(!this->isClosedShell && this->nTCS_ == 1){
      TMap EJB(this->ErrorBetaMem_ + (j%(N-1))*NBSq,NB,NB);
      TMap EKB(this->ErrorBetaMem_ + (k%(N-1))*NBSq,NB,NB);
//      B(j,k) += -EJB.frobInner(EKB);
      B(j,k) += EJB.frobInner(EKB);
    }
    B(k,j) = B(j,k);
  }

  for (auto l=0;l<N-1;l++){
     B(N-1,l)=-1.0;
     B(l,N-1)=-1.0;
  }

  B(N-1,N-1)=0;
  for(auto k = 0; k < N;k++) coef[k] = 0.0; 
  coef[N-1]=-1.0;

  TMap COEFF(coef,N,1);
  prettyPrint(this->fileio_->out,B,"CDIIS B Metric",20);
  prettyPrint(this->fileio_->out,COEFF,"CDIIS RHS");
  
  double ANORM = B.template lpNorm<1>();

  int LWORK  = 5*this->nDIISExtrap_;
  T   *WORK  = this->memManager_->template malloc<T>(LWORK);

  if(typeid(T).hash_code() == typeid(dcomplex).hash_code()){
    int LRWORK = 3*this->nDIISExtrap_;
    double   *RWORK = this->memManager_->template malloc<double>(LRWORK);
    zgesv_(&N,&NRHS,reinterpret_cast<dcomplex*>(B.data()),&N,iPiv,
        reinterpret_cast<dcomplex*>(coef),&N,&INFO);
    zgecon_(&NORM,&N,reinterpret_cast<dcomplex*>(B.data()),&N,&ANORM,&RCOND,
        reinterpret_cast<dcomplex*>(WORK),RWORK,&INFO);
    this->memManager_->free(RWORK,LRWORK);
  } else {
    dgesv_(&N,&NRHS,reinterpret_cast<double*>(B.data()),&N,iPiv,
        reinterpret_cast<double*>(coef),&N,&INFO);
    dgecon_(&NORM,&N,reinterpret_cast<double*>(B.data()),&N,&ANORM,&RCOND,
        reinterpret_cast<double*>(WORK),iWORK_,&INFO);
  }
  this->memManager_->template free(WORK,LWORK);


  if(std::abs(RCOND) > std::numeric_limits<double>::epsilon()) {
    this->fockA_->setZero();
    if(!this->isClosedShell && this->nTCS_ == 1) 
      this->fockB_->setZero();

    for(auto j = 0; j < N-1; j++) {
      TMap FA(this->FADIIS_ + (j%(N-1))*NBSq,NB,NB);
      (*this->fockA_) += coef[j]*FA;
      if(!this->isClosedShell && this->nTCS_ == 1) {
        TMap FB(this->FBDIIS_ + (j%(N-1))*NBSq,NB,NB);
        (*this->fockB_) += coef[j]*FB;
      }
    }
  } else {
    this->fileio_->out << std::setw(2) << " " <<
      std::setw(4) << " " <<
      "*** CDIIS Extrapolation Failed (RCOND = " << std::abs(RCOND) << 
      ") ***" << endl;
  }
  this->memManager_->free(B.data(),N*N);
  this->memManager_->free(coef,N);
  this->memManager_->free(iPiv,N);
  this->memManager_->free(iWORK_,N);

  if(!this->isClosedShell && this->nTCS_ !=2){
    (*this->fockScalar_) = 0.5 * ((*this->fockA_) + (*this->fockB_));
    (*this->fockMz_) = 0.5 * ((*this->fockA_) - (*this->fockB_));
  }

} // CDIIS

template<typename T>
void SingleSlater<T>::CpyFock(int iter){
  auto NTCSxNBASIS = this->nTCS_ * this->nBasis_;
  auto NSQ = NTCSxNBASIS  * NTCSxNBASIS;
  std::memcpy(this->FADIIS_+(iter % (this->nDIISExtrap_-1)) * NSQ,
    this->fockA_->data(),NSQ * sizeof(T));

  if(!this->isClosedShell && this->nTCS_ == 1)
    std::memcpy(this->FBDIIS_ + (iter % (this->nDIISExtrap_-1)) * NSQ,
                this->fockB_->data(),NSQ * sizeof(T));
} // CpyFock

template<typename T>
void SingleSlater<T>::cpyOrthoFock2(int iter){
  T* ScalarPtr;
  if(this->isClosedShell && this->nTCS_ == 1) 
    ScalarPtr = this->fockOrthoA_->data();
//  ScalarPtr = this->fockA_->data();
  else
    ScalarPtr = this->fockOrthoScalar_->data();
//  ScalarPtr = this->fockScalar_->data();

  hsize_t offset[] = {iter % (this->nDIISExtrap_-1),0,0};
  hsize_t stride[] = {1,1,1};
  hsize_t block[]  = {1,1,1};
 
  hsize_t subDim[] = {1,this->nBasis_,this->nBasis_};
  hsize_t count[]  = {1,this->nBasis_,this->nBasis_};

  cout << "COPY FOCK " << iter << endl;
  H5::DataSpace FScalar, FMz, FMy, FMx;
  FScalar = this->FScalarDIIS_->getSpace();
  FScalar.selectHyperslab(H5S_SELECT_SET,count,offset,stride,block);
  if(this->nTCS_ == 2 || !this->isClosedShell){
    FMz = this->FMzDIIS_->getSpace();
    FMz.selectHyperslab(H5S_SELECT_SET,count,offset,stride,block);
  }
  if(this->nTCS_ == 2) {
    FMy = this->FMyDIIS_->getSpace();
    FMx = this->FMxDIIS_->getSpace();
    FMy.selectHyperslab(H5S_SELECT_SET,count,offset,stride,block);
    FMx.selectHyperslab(H5S_SELECT_SET,count,offset,stride,block);
  }
  
  H5::DataSpace memSpace(3,subDim,NULL);

  this->FScalarDIIS_->write(ScalarPtr,H5PredType<T>(),memSpace,FScalar); 
  if(this->nTCS_ == 2 || !this->isClosedShell)
    this->FMzDIIS_->write(this->fockOrthoMz_->data(),H5PredType<T>(),
      memSpace,FMz); 
//  this->FMzDIIS_->write(this->fockMz_->data(),H5PredType<T>(),
//    memSpace,FMz); 
  if(this->nTCS_ == 2) {
    this->FMyDIIS_->write(this->fockOrthoMy_->data(),H5PredType<T>(),
      memSpace,FMy); 
    this->FMxDIIS_->write(this->fockOrthoMx_->data(),H5PredType<T>(),
      memSpace,FMz); 
  //this->FMyDIIS_->write(this->fockMy_->data(),H5PredType<T>(),memSpace,FMy); 
  //this->FMxDIIS_->write(this->fockMx_->data(),H5PredType<T>(),memSpace,FMz); 
  }
}

template<typename T>
void SingleSlater<T>::cpyOrthoDen2(int iter){
  T* ScalarPtr;
  if(this->isClosedShell && this->nTCS_ == 1) 
    ScalarPtr = this->onePDMOrthoA_->data();
  else
    ScalarPtr = this->onePDMOrthoScalar_->data();

  hsize_t offset[] = {iter % (this->nDIISExtrap_-1),0,0};
  hsize_t stride[] = {1,1,1};
  hsize_t block[]  = {1,1,1};
 
  hsize_t subDim[] = {1,this->nBasis_,this->nBasis_};
  hsize_t count[]  = {1,this->nBasis_,this->nBasis_};

  H5::DataSpace DScalar, DMz, DMy, DMx;
  DScalar = this->DScalarDIIS_->getSpace();
  DScalar.selectHyperslab(H5S_SELECT_SET,count,offset,stride,block);
  if(this->nTCS_ == 2 || !this->isClosedShell){
    DMz = this->DMzDIIS_->getSpace();
    DMz.selectHyperslab(H5S_SELECT_SET,count,offset,stride,block);
  }
  if(this->nTCS_ == 2) {
    DMy = this->DMyDIIS_->getSpace();
    DMx = this->DMxDIIS_->getSpace();
    DMy.selectHyperslab(H5S_SELECT_SET,count,offset,stride,block);
    DMx.selectHyperslab(H5S_SELECT_SET,count,offset,stride,block);
  }

  H5::DataSpace memSpace(3,subDim,NULL);
  
  this->DScalarDIIS_->write(ScalarPtr,H5PredType<T>(),memSpace,DScalar); 
  if(this->nTCS_ == 2 || !this->isClosedShell)
    this->DMzDIIS_->write(this->onePDMOrthoMz_->data(),H5PredType<T>(),
      memSpace,DMz); 
  if(this->nTCS_ == 2) {
    this->DMyDIIS_->write(this->onePDMOrthoMy_->data(),H5PredType<T>(),
      memSpace,DMy); 
    this->DMxDIIS_->write(this->onePDMOrthoMx_->data(),H5PredType<T>(),
      memSpace,DMx); 
  }
}

template<typename T>
void SingleSlater<T>::cpyAOtoOrthoDen(){
  if(this->nTCS_ == 1 && this->isClosedShell) {
    // Just copy for RHF
    (*this->onePDMOrthoA_) = (*this->onePDMA_);
  } else if(this->nTCS_ == 2 || !this->isClosedShell){
    // Scatter 
    std::vector<std::reference_wrapper<TMap>> scattered;
    scattered.emplace_back(*this->onePDMOrthoScalar_);
    scattered.emplace_back(*this->onePDMOrthoMz_);

    if(this->nTCS_ == 1) {
      Quantum<T>::spinScatter(*this->onePDMA_,*this->onePDMB_,scattered);
    } else {
      scattered.emplace_back(*this->onePDMOrthoMy_);
      scattered.emplace_back(*this->onePDMOrthoMx_);
      Quantum<T>::spinScatter(*this->onePDMA_,scattered);
    }

    /*
    // Copy
    (*this->onePDMOrthoScalar_) = (*this->onePDMScalar_);
    (*this->onePDMOrthoMz_)     = (*this->onePDMMz_);
    if(this->nTCS_ == 2) {
      (*this->onePDMOrthoMy_)     = (*this->onePDMMy_);
      (*this->onePDMOrthoMx_)     = (*this->onePDMMx_);
    }
    */
  }
};

template<typename F> inline F DIISComplexScale();
template<>
inline double DIISComplexScale<double>(){return -1.; }
template<>
inline dcomplex DIISComplexScale<dcomplex>(){return dcomplex(0,1); }

template<typename T>
void SingleSlater<T>::genDIISCom(int iter){
  hsize_t offset[] = {iter % (this->nDIISExtrap_-1),0,0};
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
    this->NBSqScratch_->noalias() = 
      (*this->fockOrthoA_) * (*this->onePDMOrthoA_);
    this->NBSqScratch_->noalias() -=
      (*this->onePDMOrthoA_) * (*this->fockOrthoA_);
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
