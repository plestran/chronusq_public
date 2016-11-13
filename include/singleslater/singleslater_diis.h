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

//#include <singleslater/singleslater_olddiis.h>
#include <extrapolate.h>

template<typename F> inline F DIISComplexScale();
template<>
inline double DIISComplexScale<double>(){return -1.; }
template<>
inline dcomplex DIISComplexScale<dcomplex>(){return dcomplex(0,1); }

template<typename T>
void SingleSlater<T>::CDIIS4(int NDIIS){

    std::function<void(H5::DataSet*,std::size_t,T*)> F1 = 
      std::bind(&SingleSlater<T>::readDIIS,this,
        std::placeholders::_1,std::placeholders::_2,std::placeholders::_3);

    std::function<T*(std::size_t)> F2 = 
      std::bind(&CQMemManager::malloc<T>,this->memManager_,
        std::placeholders::_1);

    std::function<void(T*,std::size_t)> F3 = 
      std::bind(&CQMemManager::free<T>,this->memManager_,
        std::placeholders::_1,std::placeholders::_2);

    T* coef;
    auto Extrap = [&](TMap *res, H5::DataSet *basis) {
      for(auto j = 0; j < NDIIS; j++) {
        this->readDIIS(basis,j,this->NBSqScratch_->data());
        (*res) += coef[j] * (*this->NBSqScratch_); 
      }
    };

    DIIS<T> extrap(NDIIS,this->nBasis_*this->nBasis_,F1,F2,F3,
      this->CommDIIS_);

    if(!extrap.extrapolate()) {
      this->fileio_->out << 
      "    *** DIIS Inversion failed: Defaulting to Conventional SCF Step ***" 
        << endl;
      return;
    };

    coef = extrap.coeffs(); 


    for(auto FOCK : this->fock_) FOCK->setZero();
    for(auto PT   : this->PT_)   PT->setZero();

    for(auto I = 0; I < this->fock_.size(); I++){
      Extrap(this->fock_[I],this->FDIIS_[I]);
      Extrap(this->PT_[I],this->PTDIIS_[I]);
    }
/*
    for(auto I = 0; I < this->fock_.size(); I++){
      std::vector<H5::DataSet*> eFiles;
      eFiles.emplace_back(this->CommDIIS_[I]);

      DIIS<T> extrap(NDIIS,this->nBasis_*this->nBasis_,F1,F2,F3,
        eFiles);

      if(!extrap.extrapolate()) continue;

      coef = extrap.coeffs(); 
      this->fock_[I]->setZero();
      this->PT_[I]->setZero();

      Extrap(this->fock_[I],this->FDIIS_[I]);
      Extrap(this->PT_[I],this->PTDIIS_[I]);
    }
*/


  this->orthoFock3();

}


template<typename T>
void SingleSlater<T>::genDIISCom(int iter){
  int ITER = iter % this->nDIISExtrap_;
  double XSmall = std::numeric_limits<double>::epsilon()*10*5;

  FPScalar_->read(this->NBSqScratch_->data(),H5PredType<T>());

  this->NBSqScratch2_->noalias() = 
    (*this->NBSqScratch_) - this->NBSqScratch_->adjoint();
  this->aointegrals_->Ortho2Trans(*this->NBSqScratch2_,*this->NBSqScratch2_);

  this->isConverged = this->NBSqScratch2_->norm() < XSmall;


  this->writeDIIS(this->EScalarDIIS_,ITER,this->NBSqScratch2_->data());

  if(this->nTCS_ == 2 or !this->isClosedShell) {
    FPMz_->read(this->NBSqScratch_->data(),H5PredType<T>());

    this->NBSqScratch2_->noalias() = 
      (*this->NBSqScratch_) - this->NBSqScratch_->adjoint();

    this->aointegrals_->Ortho2Trans(*this->NBSqScratch2_,*this->NBSqScratch2_);

    this->writeDIIS(this->EMzDIIS_,ITER,this->NBSqScratch2_->data());
    this->isConverged =  this->isConverged and
      this->NBSqScratch2_->norm() < XSmall;
  }

  if(this->nTCS_ == 2) {
    FPMy_->read(this->NBSqScratch_->data(),H5PredType<T>());

    this->NBSqScratch2_->noalias() = 
      (*this->NBSqScratch_) - this->NBSqScratch_->adjoint();
    this->aointegrals_->Ortho2Trans(*this->NBSqScratch2_,*this->NBSqScratch2_);

    this->writeDIIS(this->EMyDIIS_,ITER,this->NBSqScratch2_->data());
    this->isConverged =  this->isConverged and
      this->NBSqScratch2_->norm() < XSmall;

    FPMx_->read(this->NBSqScratch_->data(),H5PredType<T>());

    this->NBSqScratch2_->noalias() = 
      (*this->NBSqScratch_) - this->NBSqScratch_->adjoint();
    this->aointegrals_->Ortho2Trans(*this->NBSqScratch2_,*this->NBSqScratch2_);

    this->writeDIIS(this->EMxDIIS_,ITER,this->NBSqScratch2_->data());
    this->isConverged =  this->isConverged and
      this->NBSqScratch2_->norm() < XSmall;
  }

  if(this->isConverged and this->doDIIS){
    this->diisAlg_ = NO_DIIS;
  }
  
}

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

  this->FDIIS_.emplace_back(this->FScalarDIIS_);
  this->DDIIS_.emplace_back(this->DScalarDIIS_);
  this->CommDIIS_.emplace_back(this->EScalarDIIS_);
  this->PTDIIS_.emplace_back(this->PTScalarDIIS_);

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

    this->FDIIS_.emplace_back(this->FMzDIIS_);
    this->DDIIS_.emplace_back(this->DMzDIIS_);
    this->CommDIIS_.emplace_back(this->EMzDIIS_);
    this->PTDIIS_.emplace_back(this->PTMzDIIS_);
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

    this->FDIIS_.emplace_back(this->FMyDIIS_);
    this->DDIIS_.emplace_back(this->DMyDIIS_);
    this->CommDIIS_.emplace_back(this->EMyDIIS_);
    this->PTDIIS_.emplace_back(this->PTMyDIIS_);

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

    this->FDIIS_.emplace_back(this->FMxDIIS_);
    this->DDIIS_.emplace_back(this->DMxDIIS_);
    this->CommDIIS_.emplace_back(this->EMxDIIS_);
    this->PTDIIS_.emplace_back(this->PTMxDIIS_);
  }

};


template<typename T>
void SingleSlater<T>::cpyDenDIIS(int iter){
  int ITER = iter % this->nDIISExtrap_;
  for(auto I = 0; I < this->onePDM_.size(); I++){
    this->writeDIIS(this->DDIIS_[I],ITER,this->onePDMOrtho_[I]->data());
  }
};

template<typename T>
void SingleSlater<T>::cpyFockDIIS(int iter){
  int ITER = iter % this->nDIISExtrap_;
  for(auto I = 0; I < this->fock_.size(); I++){
    this->writeDIIS(this->FDIIS_[I],ITER,this->fock_[I]->data());
    this->writeDIIS(this->PTDIIS_[I],ITER,this->PT_[I]->data());
  }
};

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
